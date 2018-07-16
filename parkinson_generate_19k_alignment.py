import pandas as pd
import os
from multiprocessing import Queue, Process, Manager, current_process
from plumbum import local
import os
import sys
import subprocess
from collections import defaultdict
def readDefinedFileToList(filename):
    temp_list = []
    with open(filename, mode='r') as reader:
        temp_list = [line.rstrip() for line in reader]
    return temp_list

def writeListToDestination(destination, listToWrite):
    #print('Writing list to ' + destination)
    try:
        os.makedirs(os.path.dirname(destination))
    except FileExistsError:
        pass

    with open(destination, mode='w') as writer:
        i = 0
        while i < len(listToWrite):
            if i != len(listToWrite)-1:
                writer.write(listToWrite[i] + '\n')
            elif i == len(listToWrite)-1:
                writer.write(listToWrite[i])
            i += 1



def convert_interleaved_to_sequencial_fasta_two(fasta_in):
    fasta_out = []

    for i in range(len(fasta_in)):

        if fasta_in[i].startswith('>'):
            if fasta_out:
                # if the fasta is not empty then this is not the first
                fasta_out.append(temp_seq_str)
            #else then this is the first sequence and there is no need to add the seq.
            temp_seq_str = ''
            fasta_out.append(fasta_in[i])
        else:
            temp_seq_str = temp_seq_str + fasta_in[i]
    #finally we need to add in the last sequence
    fasta_out.append(temp_seq_str)
    return fasta_out


def generate_local_alignments_for_each_ortholog():
    # the amino acid sequences
    aa_seq_array = pd.read_csv('/home/humebc/projects/parky/aa_seq_multi_orf_orths_fixed.csv', sep=',', lineterminator='\n', index_col=0, header=0)

    # the gene IDs
    gene_id_array = pd.read_csv('/home/humebc/projects/parky/gene_id_fixed.csv', sep=',', lineterminator='\n', index_col=0, header=0)

    # the row and column headings should be the same for both arrays so we can iterate both at the same time
    # for each row
    # each row represents an ortholog. These will be our units to process. We may want to multiprocess this.
    # to make this multiprocess friendly, lets throw tuples of gene_id_name and aa_sequence into a list and then
    # add this to a list. Creating a 2D list. We can then offer up each of these list items to a queue that we can parse
    # when doing the MPing

    # making MP data_holder_list
    tuple_holder_dict = {}
    #for each ortholog
    for row_index in aa_seq_array.index.values.tolist():
        print('Adding ortholog {} to MP info list'.format(row_index))
        ortholog_id_seq_list = []
        # for each species.
        for spp in list(gene_id_array):
            gene_id = gene_id_array[spp][row_index]
            aa_seq = aa_seq_array[spp][row_index]
            ortholog_id_seq_list.append((gene_id, aa_seq))
        tuple_holder_dict[row_index] = ortholog_id_seq_list

    # creating the MP input queue
    ortholog_input_queue = Queue()

    # populate with one key value pair per ortholog
    for key, value in tuple_holder_dict.items():
        print('Placing {} in MP queue'.format(key))
        ortholog_input_queue.put((key, value))

    num_proc = 24

    # put in the STOPs
    for N in range(num_proc):
        ortholog_input_queue.put('STOP')

    allProcesses = []

    # directory to put the local alignments
    output_dir = '/home/humebc/projects/parky/aa_tree_creation/local_alignments'

    # the list of species for each ortholog
    spp_list = [spp for spp in list(gene_id_array)]

    # Then start the workers
    for N in range(num_proc):
        p = Process(target=create_local_alignment_worker, args=(ortholog_input_queue, output_dir, spp_list))
        allProcesses.append(p)
        p.start()

    for p in allProcesses:
        p.join()

    # at this point we have the local alignments all written as fasta files to output_dir.
    # Now it just remains to concatenate and then run the ML tree.


def create_local_alignment_worker(input, output_dir, spp_list):
    # for each list that represents an ortholog
    for k_v_pair in iter(input.get, 'STOP'):

        # ortholog_id
        ortholog_id = k_v_pair[0]
        print('Processing {}'.format(ortholog_id))
        # ortholog_spp_seq_list
        ortholog_spp_seq_list = k_v_pair[1]

        # create the fasta
        fasta_file = []

        # for each species
        # add the name and aa_seq to the fasta_file
        for i in range(len(spp_list)):
            fasta_file.extend(['>{}_{}'.format(spp_list[i], ortholog_spp_seq_list[i][0]), ortholog_spp_seq_list[i][1]])

        # here we have the fasta_file populated

        # Write out the new fasta
        fasta_output_path = '{}/{}.fasta'.format(output_dir, ortholog_id)
        writeListToDestination(fasta_output_path, fasta_file)
        # now perform the alignment with MAFFT
        mafft = local["mafft-linsi"]
        in_file = fasta_output_path
        out_file = fasta_output_path.replace('.fasta', '_aligned.fasta')
        # now run mafft including the redirect
        (mafft[in_file] > out_file)()

        # at this point we have the aligned .fasta written to the output directory
        # at this point we need to trim the fasta.
        # I was going to use trimAl but this doesn't actually have an option to clean up the ends of alignments
        # instead, read in the alignment as a TwoD list to a pandas dataframe
        # then delete the begining and end columns that contain gap sites
        aligned_fasta_interleaved = readDefinedFileToList(out_file)
        aligned_fasta = convert_interleaved_to_sequencial_fasta_two(aligned_fasta_interleaved)
        array_list = []
        for i in range(1, len(aligned_fasta), 2):
                array_list.append(list(aligned_fasta[i]))

        # make into pandas dataframe
        alignment_df = pd.DataFrame(array_list)

        # go from either end deleting any columns that have a gap
        columns_to_drop = []
        for i in list(alignment_df):
            # if there is a gap in the column at the beginning
            if '-' in list(alignment_df[i]) or '*' in list(alignment_df[i]):
                columns_to_drop.append(i)
            else:
                break
        for i in reversed(list(alignment_df)):
            # if there is a gap in the column at the end
            if '-' in list(alignment_df[i]) or '*' in list(alignment_df[i]):
                columns_to_drop.append(i)
            else:
                break

        # get a list that is the columns indices that we want to keep
        col_to_keep = [col_index for col_index in list(alignment_df) if col_index not in columns_to_drop]

        # drop the gap columns
        alignment_df = alignment_df[col_to_keep]

        # here we have the pandas dataframe with the gap columns removed
        #convert back to a fasta and write out
        cropped_fasta = []
        alignment_index_labels = list(alignment_df.index)
        for i in range(len(alignment_index_labels)):
            seq_name = '>{}_{}'.format(spp_list[i], ortholog_spp_seq_list[i][0])
            aa_seq = ''.join(alignment_df.loc[alignment_index_labels[i]])
            cropped_fasta.extend([seq_name, aa_seq])

        # here we have the cropped and aligned fasta
        # write it out
        aligned_cropped_fasta_path = fasta_output_path.replace('.fasta', '_aligned_cropped.fasta')
        writeListToDestination(aligned_cropped_fasta_path, cropped_fasta)

        # here we should be done with the single alignment
        print('Local alignment for {} completed'.format(ortholog_id))

def concatenate_local_alignments():
    # The master alignment that we create should be partitioned according to the protein model used.
    # I have generated all of the .out files which are the outputs from the prottest
    # We should iter through these and create a dictionary that is a model type as key
    # and then have a list of the orthologs of that model.
    # then sort this by the length of the list
    # then work our way through the local alignments in this order creating the alignment
    # We will need to generate the p file as we go
    # this should take the form
    '''
    JTT, gene1 = 1-500
    WAGF, gene2 = 501-800
    WAG, gene3 = 801-1000

    '''

    # get list of the .out prottest files
    base_dir = '/home/humebc/projects/parky/aa_tree_creation/local_alignments'
    list_of_prot_out_filenames = [f for f in os.listdir(base_dir) if 'prottest' in f]

    # iter through the list of protfiles creating the dict relating model to ortholog
    # we cannnot change the +G or +I for each partition. As such I will define according to the base model
    model_to_orth_dict = defaultdict(list)
    for i in range(len(list_of_prot_out_filenames)):
        model = ''
        file_name = list_of_prot_out_filenames[i]
        orth_num = int(file_name.split('_')[0])
        with open('{}/{}'.format(base_dir, file_name), 'r') as f:
            temp_file_list = [line.rstrip() for line in f]
        for j in range(300, len(temp_file_list), 1):
            if 'Best model according to BIC' in temp_file_list[j]:
                model = temp_file_list[j].split(':')[1].strip().replace('+G','').replace('+I','')
                break
        if model == '':
            sys.exit('Model line not found in {}'.format(orth_num))
        model_to_orth_dict[model].append(orth_num)

    # #N.B. that we cannot have different gamma for different partitions
    # # same for +I so I will run with GAMMAI as the -m argument
    # for model in model_to_orth_dict

    print('The 19k sequences are best represented by {} different aa models'.format(len(model_to_orth_dict.keys())))

    # here we have the dict populated
    # now sort the dict
    sorted_model_list = sorted(model_to_orth_dict, key=lambda k: len(model_to_orth_dict[k]), reverse=True)

    # now go model by model in the sorted_model_list to make the master alignment.

    # get a list of all of the fasta names that we will want to concatenate

    # list_of_fasta_file_names = [f for f in os.listdir(base_dir) if 'aligned_cropped.fasta' in f]

    # # lets sort the list of files by number so that we start with the smallest
    # nums_in_file_names = sorted([int(name.split('_')[0]) for name in list_of_fasta_file_names])

    # list_of_files = [str(file_num) + '_aligned_cropped.fasta' for file_num in nums_in_file_names]

    # not the most elegant way but I think I'll just create the mast fasta in memory
    master_fasta = ['>min','', '>pmin', '', '>psyg', '', '>ppsyg', '']

    # The q file will hold the information for the partitioning of the alignment for the raxml analysis
    q_file = []
    for model in sorted_model_list:
        q_file_start = len(master_fasta[1]) + 1
        sys.stdout.write('\rProcessing model {} sequences'.format(model))
        for orth_num in model_to_orth_dict[model]:
            file_name = str(orth_num) + '_aligned_cropped.fasta'
            with open('{}/{}'.format(base_dir, file_name), 'r') as f:
                temp_list_of_lines = [line.rstrip() for line in f]

            for i in range(1, len(temp_list_of_lines), 2):
                new_seq = master_fasta[i] + temp_list_of_lines[i]
                master_fasta[i] = new_seq
        q_file_finish = len(master_fasta[1])
        q_file.append('{}, gene{} = {}-{}'.format(
            model.upper(), sorted_model_list.index(model) + 1, q_file_start, q_file_finish))

    # here we have the master fasta and the q file ready to be outputted

    # now write out the master fasta
    master_fasta_output_path = '/home/humebc/projects/parky/aa_tree_creation/master.fasta'
    with open(master_fasta_output_path, 'w') as f:
        for line in master_fasta:
            f.write('{}\n'.format(line))

    # now write out the q file
    q_file_output_path = '/home/humebc/projects/parky/aa_tree_creation/qfile.q'
    with open(q_file_output_path, 'w') as f:
        for line in q_file:
            f.write('{}\n'.format(line))

    # now run raxml
    #NB note that although we are specificing mdels for each partition, we still need to add the PROTGAMMAIWAG
    # model argument to the -m flag. This is just a weird operation of raxml and is explained but hidden in the manual
    # (search for 'If you want to do a partitioned analysis of concatenated'). Raxml will only extract the rate
    # variation information from this and will ignore the model component e.g. WAG. FYI any model could be used
    # doesn't have to be WAG.
    raxml_path = '/home/humebc/phylogeneticsSoftware/raxml/standard-RAxML/raxmlHPC-PTHREADS-AVX'
    subprocess.run([raxml_path, '-s', master_fasta_output_path, '-q', q_file_output_path,
                    '-x', '183746', '-f', 'a', '-p', '83746273', '-#', '100', '-T', '8', '-n', 'parkinson_out',
                    '-m', 'PROTGAMMAWAG', '-w', '/home/humebc/projects/parky/aa_tree_creation'])

    print('\nConstruction of master fasta complete:\n{}\n{}'.format(master_fasta_output_path, q_file_output_path))

def convert_fasta_to_phylip():
    from Bio import AlignIO

    base_dir = '/home/humebc/projects/parky/aa_tree_creation/local_alignments'
    fasta_alignment = AlignIO.read('{}/master.fasta'.format(base_dir), 'fasta')
    AlignIO.write(fasta_alignment, handle='/home/humebc/projects/parky/aa_tree_creation/master.phylip', format='phylip')

def generate_protein_substitution_models_for_each_gene():
    # we will find each of the local alignments and run put them into a list which we will MP
    # for each item we will run prottest with a single thread and output a file
    # in the concatenate local alignments file we will then create a q file that will
    # designate the different partitions and which of the substitution models to use.

    # get a list of all of the fasta names that we will want to concatenate
    base_dir = '/home/humebc/projects/parky/aa_tree_creation/local_alignments'
    list_of_files = [f for f in os.listdir(base_dir) if 'aligned_cropped.fasta' in f]


    num_proc = 12

    # Queue that will hold the index of the rows that need to be checked
    input_queue = Queue()

    # populate input_queue
    for file_name in list_of_files:
        input_queue.put(file_name)

    for i in range(num_proc):
        input_queue.put('STOP')

    list_of_processes = []
    for N in range(num_proc):
        p = Process(target=prottest_worker, args=(input_queue,))
        list_of_processes.append(p)
        p.start()

    for p in list_of_processes:
        p.join()

    return

def prottest_worker(input_queue):
    base_dir = '/home/humebc/projects/parky/aa_tree_creation/local_alignments'
    for file_name in iter(input_queue.get, 'STOP'):
        input_path = '{}/{}'.format(base_dir, file_name)
        output_path = input_path.replace('_aligned_cropped.fasta', '_prottest_result.out')
        if os.path.isfile(output_path):
            continue
        sys.stdout.write('\rRunning prottest: {}'.format(file_name))
        # perform the prottest
        prot_path = '/home/humebc/phylogeneticsSoftware/protest/prottest-3.4.2/prottest-3.4.2.jar'
        subprocess.run(['java', '-jar', prot_path, '-i', input_path, '-o', output_path, '-all-distributions', '-all']
                       , stdout=subprocess.PIPE, stderr=subprocess.PIPE)



concatenate_local_alignments()
























