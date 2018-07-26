import pandas as pd
import os
from multiprocessing import Queue, Process, Manager, current_process
from plumbum import local
import os
import sys
import subprocess
from collections import defaultdict
import itertools
import pickle
import shutil
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

def generate_local_alignments_for_each_ortholog_two():
    # So we were generating an alignment for the aa seqs using mafft and then cropping them
    # However, if we are going to use neighbour to do the aligning then that automatically produces us
    # a MAFFT alignment and allows us to screen for low scoring site that can additionally be translated
    # into removal of uncertain columns.
    # similar to how we were creating a directory for each of the aa_seqs we will do the same thing for
    # each of the orthologs.

    # the cds sequences
    cds_seq_df = pd.read_csv('/home/humebc/projects/parky/cds_seq_multi_orf_orths__partial_cdsfixed.csv', sep=',', lineterminator='\n',
                             index_col=0, header=0)

    # the list of species for each ortholog
    spp_list = [spp for spp in list(cds_seq_df)]


    #we will need a directory for each of the orthologs
    # in each directory we will have a cds fasta
    # to MP this we can simply give a list of the directories and get the ortholog ID from the directory name

    # for each ortholog, create directory and write in the fasta of the cds then pass the directory to the MP list

    MP_list = []


    for row_index in cds_seq_df.index.values.tolist():
        temp_fasta = []
        sys.stdout.write('\rGenerating fasta for ortholog {}'.format(row_index))

        # populate the fasta
        list_of_seq_names = []
        for spp in spp_list:
            temp_fasta.extend(['>{}_{}'.format(row_index, spp), cds_seq_df.loc[row_index, spp]])
            list_of_seq_names.append('{}_{}'.format(row_index, spp))
        # write out the fasta into the  directory
        base_dir = '/home/humebc/projects/parky/local_alignments/{0}/'.format(row_index)
        path_to_fasta = '/home/humebc/projects/parky/local_alignments/{0}/unaligned_cds_{0}.fasta'.format(row_index)
        os.makedirs(base_dir, exist_ok=True)
        with open(path_to_fasta, 'w') as f:
            for line in temp_fasta:
                f.write('{}\n'.format(line))

        # add the fasta full path to the MP list
        MP_list.append((path_to_fasta, list_of_seq_names))




    # creating the MP input queue
    ortholog_input_queue = Queue()

    # populate with one key value pair per ortholog
    for fasta_tup in MP_list:
        ortholog_input_queue.put(fasta_tup)

    num_proc = 24

    # put in the STOPs
    for N in range(num_proc):
        ortholog_input_queue.put('STOP')

    allProcesses = []


    # Then start the workers
    for N in range(num_proc):
        p = Process(target=guidance_worker, args=(ortholog_input_queue,))
        allProcesses.append(p)
        p.start()

    for p in allProcesses:
        p.join()

        # here we have the cds and the aa alignments cropped and written out
        # we can then use these as input into CODEML and the BUSTED programs
        # and to run the model calls for making the tree with.

def guidance_worker(input_queue):
    for fasta_tup in iter(input_queue.get, 'STOP'):
        fasta_file_full_path, seq_names = fasta_tup
        ortholog_id = fasta_file_full_path.split('/')[-2]
        sys.stdout.write('\rPerforming Guidance for {}'.format(ortholog_id))
        guidance_full_path = '/home/humebc/phylogeneticsSoftware/guidance2/guidance.v2.02/www/Guidance/guidance.pl'
        output_dir_full_path = '/'.join(fasta_file_full_path.split('/')[:-1])

        # check to see if guidance has already been completed

        # if os.path.isfile('{}/{}.cropped_aligned_cds.fasta'.format(output_dir_full_path, ortholog_id)):
        #     continue

        # subprocess.run(['perl', guidance_full_path, '--seqFile',
        #                 fasta_file_full_path, '--msaProgram', 'MAFFT', '--seqType',
        #                 'codon', '--outDir', output_dir_full_path, '--bootstraps', '10', '--outOrder', 'as_input',
        #                 '--colCutoff', '0.6', '--dataset', ortholog_id], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # # at this point we should have performed the guidance analysis
        # # we can now read in the cols that should be removed from the cds and aa alignments
        # cols_to_remove_cds = []
        # cds_cols_to_remove_file_path = 'Seqs.Orig_DNA.fas.FIXED.{}.MAFFT.Removed_Col'.format(ortholog_id)
        # with open(cds_cols_to_remove_file_path, 'r') as f:
        #     cds_col_file_list = [line.rstrip() for line in f]
        # for line in cds_col_file_list:
        #     cols_to_remove_cds.append(int(line.split()[2]))

        # at this point we should have performed the guidance analysis
        # we can now read in the scores for each aa residue and drop any alignment columns
        # that have only -na scores or contain a score < the cutoff which is 0.6
        # we will be working in sets of four because there are four sequences. Each column
        # we find to delte in this aa score matrix will represent three columns to be dropped
        # in the cds alignment.
        # also, bear in mind that the columns of the alignment are not 0 indexed but start at 1
        # we will have to take this into account when working out which columns of the alignments to drop

        cols_to_remove_aa = []
        aa_cols_score_file_path = '{}/{}.MAFFT.Guidance2_res_pair_res.PROT.scr'.format(output_dir_full_path, ortholog_id)

        with open(aa_cols_score_file_path, 'r') as f:
            aa_col_score_file_list = [line.rstrip() for line in f][1:]

        for i in range(int(len(aa_col_score_file_list)/4)):

            socres_set = []
            # get a list of the scores for each position
            # for j in range(i, i+4, 1):
            #     scores_set.append()
            # the range that represents i and the next three i s in the iteration without increasing i
            file_line_indices = [i*4 + k for k in range(4)]
            scores_set = [float(aa_col_score_file_list[j].split()[2]) if '-nan' not in aa_col_score_file_list[j] else '-nan' for j in file_line_indices]

            # now examine what is in the score sets

            drop = False
            nan_count = 0
            for score in scores_set:
                if score != '-nan':
                    if score < 0.6:
                        drop = True
                        break
                elif score == '-nan':
                    nan_count += 1
                else:
                    continue
            # when we come out ned to check if the nan_score == 4
            if nan_count == 4:
                drop = True


            # if drop = True then this is a column to drop
            if drop:
                cols_to_remove_aa.append(i)

        # here we have a 0 based list of the columns that need to be dropped from the aa alignment
        # convert this to a 0 based index of the columns that need to be dropped from the cds alignment
        cols_to_remove_cds = []
        for col_index in cols_to_remove_aa:
            cols_to_remove_cds.extend([col_index*3 + n for n in range(3)])

        # here we have the indices that need to be dropped for the cds and aa alignments
        # now read in the alignments as 2D lists, convert to pandas dataframe and perform the columns drops

        # aa first
        aa_alignment_file_path = '{}/{}.MAFFT.PROT.aln'.format(output_dir_full_path, ortholog_id)
        with open(aa_alignment_file_path, 'r') as f:
            aa_alignment_file_list = convert_interleaved_to_sequencial_fasta_two([line.rstrip() for line in f])
        aa_df = pd.DataFrame([list(line) for line in aa_alignment_file_list if not line.startswith('>')])

        # now drop the columns
        columns_to_keep_aa = [col for col in list(aa_df) if col not in cols_to_remove_aa]
        aa_df = aa_df[columns_to_keep_aa]

        # cds second
        cds_alignment_file_path = '{}/{}.MAFFT.aln'.format(output_dir_full_path, ortholog_id)
        with open(cds_alignment_file_path, 'r') as f:
            cds_alignment_file_list = convert_interleaved_to_sequencial_fasta_two([line.rstrip() for line in f])
        cds_df = pd.DataFrame([list(line) for line in cds_alignment_file_list if not line.startswith('>')])

        # now drop the columns
        columns_to_keep_cds = [col for col in list(cds_df) if col not in cols_to_remove_cds]
        cds_df = cds_df[columns_to_keep_cds]

        # here we have the cds and aa dfs that we can now do the cropping with and then finally write back out as
        # fasta files
        # go from either end deleting any columns that have a gap

        # aa first
        aa_df = crop_fasta_df(aligned_fasta_as_pandas_df_to_crop = aa_df)

        # cds second
        cds_df = crop_fasta_df(aligned_fasta_as_pandas_df_to_crop = cds_df)

        # now we just need to write out dfs
        aa_fasta_list = pandas_df_to_fasta(pd_df=aa_df, seq_names=seq_names)
        with open('{}/{}.cropped_aligned_aa.fasta'.format(output_dir_full_path, ortholog_id), 'w') as f:
            for line in aa_fasta_list:
                f.write('{}\n'.format(line))

        cds_fasta_list = pandas_df_to_fasta(pd_df=cds_df, seq_names=seq_names)
        with open('{}/{}.cropped_aligned_cds.fasta'.format(output_dir_full_path, ortholog_id), 'w') as f:
            for line in cds_fasta_list:
                f.write('{}\n'.format(line))

        # here we have the cds and the aa alignments cropped and written out
        # we can then use these as input into CODEML and the BUSTED programs


def pandas_df_to_fasta(pd_df, seq_names):
    temp_fasta = []
    for i in range(len(seq_names)):
        temp_fasta.extend(['>{}'.format(seq_names[i]), ''.join(list(pd_df.loc[i]))])
    return temp_fasta


def crop_fasta_df(aligned_fasta_as_pandas_df_to_crop):
    columns_to_drop = []
    for i in list(aligned_fasta_as_pandas_df_to_crop):
        # if there is a gap in the column at the beginning
        if '-' in list(aligned_fasta_as_pandas_df_to_crop[i]) or '*' in list(aligned_fasta_as_pandas_df_to_crop[i]):
            columns_to_drop.append(i)
        else:
            break
    for i in reversed(list(aligned_fasta_as_pandas_df_to_crop)):
        # if there is a gap in the column at the end
        if '-' in list(aligned_fasta_as_pandas_df_to_crop[i]) or '*' in list(aligned_fasta_as_pandas_df_to_crop[i]):
            columns_to_drop.append(i)
        else:
            break

    # get a list that is the columns indices that we want to keep
    col_to_keep = [col_index for col_index in list(aligned_fasta_as_pandas_df_to_crop) if col_index not in columns_to_drop]
    # drop the gap columns
    return aligned_fasta_as_pandas_df_to_crop[col_to_keep]


def generate_local_alignments_for_each_ortholog():
    #TODO So we were generating an alignment for the aa seqs using mafft and then cropping them
    #However, if we are going to use neighbour to do the aligning then that automatically produces us
    # a MAFFT alignment and allows us to screen for low scoring site that can additionally be translated
    # into removal of uncertain columns.
    # similar to how we were creating a directory for each of the aa_seqs we will do the same thing for
    # each of the orthologs.


    # the amino acid sequences
    aa_seq_array = pd.read_csv('/home/humebc/projects/parky/aa_seq_multi_orf_orths_fixed.csv', sep=',', lineterminator='\n', index_col=0, header=0)

    # the gene IDs
    gene_id_array = pd.read_csv('/home/humebc/projects/parky/gene_id_fixed.csv', sep=',', lineterminator='\n', index_col=0, header=0)

    # the cds sequences
    cds_seq_df = pd.read_csv('/home/humebc/projects/parky/gene_id_fixed.csv', sep=',', lineterminator='\n',
                                index_col=0, header=0)

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
    # # Also best advice is not to run +G and +I together.
    # # As such we only need to extract the base model here i.e. WAG rather than WAG [+G|+I]
    # for model in model_to_orth_dict

    print('The 19k sequences are best represented by {} different aa models'.format(len(model_to_orth_dict.keys())))

    # here we have the dict populated
    # now sort the dict
    sorted_model_list = sorted(model_to_orth_dict, key=lambda k: len(model_to_orth_dict[k]), reverse=True)

    # now go model by model in the sorted_model_list to make the master alignment.

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
                    '-x', '183746', '-f', 'a', '-p', '83746273', '-#', '1000', '-T', '8', '-n', 'parkinson_out',
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



def generate_block_phylip_alignments_for_CODEML():
    ''' The documentation for what format the files should be in for submitting to PAML/CODEML are not so
    great. But from some testing and looking through the examples that are available two things seem to be key
    1 - automatic pairwise comparisons of sequences can be performed using the runmode as -2
    2 - a blocked alignment format is the easiest way to test all of the orthologs at once but unlike
    in the documentation, (using the G marker) this is easiest done by simply having a phymil alignment format
    one after another in succestion in a sinlge file and then setting the ndata setting to how ever many
    sets of aligments you have (about 19K for us). I've run some tests and this all seems to be running well.
    Unfortunately, the only way to MP this seems to be physically splitting up the data and running multiple
    instance of the CODEML executable.
    So... with this in mind, I will split up the 19K or so sequences into phylip alignments of 1000 sequences.
    I will then make respective control files for each of these, put them in their own directory and then set
    an instance of the CODEML running for each of these.'''


    spp_list = ['min', 'pmin', 'psyg', 'ppsyg']
    # for each pairwise comparison of the species
    wkd = '/home/humebc/projects/parky/local_alignments'
    output_dir = '/home/humebc/projects/parky/guidance_analyses'
    list_of_guidance_dirs = []
    # set that will hold the ID of any orthologs that have dud alignments
    len_zero_dir_list = []
    # counter that we'll use to split into 20 roughly 1k sequence alignments
    counter = 0
    block_counter = 0

    # for each directory or ortholog
    list_of_dirs = list()
    for root, dirs, files in os.walk(wkd):
        list_of_dirs = dirs
        break
    for dir in list_of_dirs:

        sys.stdout.write('\rOrtholog: {}     {}/{}'.format(dir, counter, len(list_of_dirs)))

        # make a new alignment file for every 1000 individual alignments
        if counter % 1000 == 0:
            if block_counter != 0:

                # then we already have a block that needs writing
                seq_file_ctrl_file_tup = write_out_cntrl_and_seq_file(block_counter=block_counter, output_dir=output_dir,
                                             phylip_alignment=phylip_alignment, num_align=1000)
                list_of_guidance_dirs.append(seq_file_ctrl_file_tup)
            # once the old block is written out start a new one
            block_counter += 1
            os.makedirs('{}/block_{}'.format(output_dir, block_counter), exist_ok=True)
            print('\n\nStarting block {}'.format(block_counter))
            phylip_alignment = []

        # add single to block master alignment
        # if the fasta was empty then this will return False
        single_phylip = generate_phylip_from_fasta(dir, wkd)

        if single_phylip:
            phylip_alignment.extend(single_phylip)
            counter += 1
        else:
            # if the fasta was empty then log this and don't add anything to the counter
            len_zero_dir_list.append(dir)
    # write out a list of the poorl alignement orfs
    with open('post_guidance_0_len_cds_alignments_orthologs.txt', 'w') as f:
        for line in len_zero_dir_list:
            f.write('{}\n'.format(line))

    # now write out the final block of alignments
    seq_file_ctrl_file_tup = write_out_cntrl_and_seq_file(block_counter, output_dir, phylip_alignment, num_align=len(list_of_dirs)-len(len_zero_dir_list)-(1000*(block_counter-1)))
    list_of_guidance_dirs.append(seq_file_ctrl_file_tup)

    pickle.dump(list_of_guidance_dirs, open( '{}/list_of_guidance_dirs.pickle'.format(output_dir), "wb" ))

def generate_phylip_from_fasta(dir, wkd):
    temp_str = str()
    temp_list = list()

    with open('{}/{}/{}.cropped_aligned_cds.fasta'.format(wkd, dir, dir), 'r') as f:
        fasta_file = [line.rstrip() for line in f]

    if len(fasta_file[1]) == 0:
        # if the fasta is empty then we need to log this outside
        return False
    else:
        for line in fasta_file:
            if line.startswith('>'):
                temp_str = line[1:]
            else:
                temp_list.append('{}    {}'.format(temp_str, line))
        # finally put in the header file
        temp_list.insert(0, '\t{} {}'.format(len(temp_list), len(fasta_file[1])))

        return temp_list


def write_out_cntrl_and_seq_file(block_counter, output_dir, phylip_alignment, num_align):
    # write out the control file specific to this alignment
    ctrl_file_path = write_out_control_file(
        output_dir='{}/block_{}'.format(output_dir, block_counter),
        num_alignments=num_align,
        grp=block_counter)
    # write out the phylip file
    seq_file_path = '{}/block_{}/block_{}_cds.phylip'.format(output_dir, block_counter, block_counter)
    with open(seq_file_path, 'w') as f:
        for line in phylip_alignment:
            f.write('{}\n'.format(line))
    # write out the tree file
    tree_file = '(ppsyg:0.01524804457090833502,(min:0.00305561548082329418,pmin:0.00350296114601793013)' \
                ':0.03350192310501232812,psyg:0.01618135662493049715);'
    tree_file_path = '{}/block_{}/block_{}_tree.nwk'.format(output_dir, block_counter, block_counter)
    with open(tree_file_path, 'w') as f:
        f.write('{}\n'.format(tree_file))

    return (seq_file_path, ctrl_file_path, tree_file_path)

def write_out_control_file(output_dir, num_alignments, grp):
    seq_file_path = '{}/block_{}_cds.phylip'.format(output_dir, grp)
    out_file_path = '{}/block_{}_guidance_results.out'.format(output_dir, grp)
    ctrl_path     = '{}/block_{}_cds.ctrl'.format(output_dir, grp)
    tree_file_path = '{}/block_{}/block_{}_tree.nwk'.format(output_dir, grp, grp)
    ctrl_file = [
    'seqfile = {}'.format(seq_file_path),
    'treefile = {}'.format(tree_file_path),
    'outfile = {}'.format(out_file_path),
    'runmode = -2',
    'seqtype = 1',
    'CodonFreq = 2',
    'ndata = {}'.format(num_alignments),
    'clock = 0',
    'model = 0',
    'NSsites = 0',
    'icode = 0',
    'fix_omega = 0',
    'omega = .4',
    'cleandata = 0'
    ]

    with open(ctrl_path, 'w') as f:
        for line in ctrl_file:
            f.write('{}\n'.format(line))

    return ctrl_path

def run_CODEML_analyses():
    ''' We should read in the tuples that contain the seq files and cntrl file
    from the pickled files that were written out from generate_master_phlip_alignments_for_CODEML
    We should then start an MP list with each of these tuples in
    In the worker we should go to each directory and start an anlysis
    The 20000 sequences were broken into 20 chunks so we should start 20 subprocess instances and
    have a one CODEML analysis run in each'''

    tup_of_dirs_list = pickle.load( open( '/home/humebc/projects/parky/guidance_analyses/list_of_guidance_dirs.pickle', "rb" ))

    input_queue = Queue()

    for tup in tup_of_dirs_list:
        input_queue.put(tup)

    num_proc = 1

    for i in range(num_proc):
        input_queue.put('STOP')

    list_of_processes = []
    for N in range(num_proc):
        p = Process(target=CODEML_run_worker, args=(input_queue,))
        list_of_processes.append(p)
        p.start()

    for p in list_of_processes:
        p.join()


def CODEML_run_worker(input_queue):
    CODEML_path = '/home/humebc/phylogeneticsSoftware/paml/paml4.9h/bin/codeml'
    for dir_tup in iter(input_queue.get, 'STOP'):
        seq_file_path, ctrl_file_path, tree_file_path = dir_tup

        wkd = os.path.dirname(seq_file_path)

        os.chdir(wkd)

        subprocess.run([CODEML_path, ctrl_file_path])


def create_directories_and_run_busted():
    ''' The BUSTED annoyingly has an interactive interface. It does accept special batch files, but I really
    don't want to have to invest the time in learning how to use those just for this analysis.
    Instead I will simply use a text file that had the predetermined answers to the interactive prompts in
    order to start analyses.
    To conduct the analyses we will have to go ortholog by ortholog, pulling out the fasta alignment
    and then writing this to a new directory were we will also write the tree. Then we simply write
    the answer file and run the command.

    We should definitely MP this. To do this we should read in the directories
    and put these in as the MP list arguments. We should then make all of the files etc and run the analyses
    within these threads.

    we need to remember to check that alignments aren't empty.'''

    local_dirs_dir = '/home/humebc/projects/parky/local_alignments'
    list_of_dirs = list()
    for root, dirs, files in os.walk(local_dirs_dir):
        list_of_dirs = dirs
        break

    input_queue = Queue()

    for dir in list_of_dirs:
        input_queue.put(dir)

    num_proc = 40

    for i in range(num_proc):
        input_queue.put('STOP')

    list_of_processes = []
    for N in range(num_proc):
        p = Process(target=BUSTED_MP_worker, args=(input_queue,))
        list_of_processes.append(p)
        p.start()

    for p in list_of_processes:
        p.join()

    apples = 'asdf'

def BUSTED_MP_worker(input_queue):
    busted_analyses_dir = '/home/humebc/projects/parky/busted_analyses'
    local_dirs_dir = '/home/humebc/projects/parky/local_alignments'


    # The HYPHYMP executable relies on some sort of relative file but I can't find out which one
    # as such we will have to change dir to the HYPHY_idr in order to be able to invoke the program.
    HYPHYMP_path = '/home/humebc/phylogeneticsSoftware/hyphy/hyphy-2.3.13/HYPHYMP'
    HYPHY_dir = '/home/humebc/phylogeneticsSoftware/hyphy/hyphy-2.3.13'
    for ortholog_id in iter(input_queue.get, 'STOP'):

        # First check to see if the busted processing has already been done for this ortholog
        # if it has then skip onto the next one.
        if os.path.isfile('{0}/{1}/{1}_qced_cds_aligned.fasta.BUSTED.json'.format(busted_analyses_dir, ortholog_id)):
            continue

        # First check to see that the alignment is good
        orig_align_path = '{0}/{1}/{1}.cropped_aligned_cds.fasta'.format(local_dirs_dir, ortholog_id)
        with open(orig_align_path, 'r') as f:
            orig_fasta = [line.rstrip() for line in f]

        if len(orig_fasta[1]) < 2:
            # then the fasta is bad and we should move onto the next directory and ignore this one.
            continue


        sys.stdout.write('\rRunning BUSTED analysis on ortholog: {}'.format(ortholog_id))
        tree_file = '({0}_ppsyg:0.01524804457090833502,({0}_min:0.00305561548082329418,{0}_pmin:0.00350296114601793013)' \
                    ':0.03350192310501232812,{0}_psyg:0.01618135662493049715);'.format(ortholog_id)

        wkd = '{}/{}'.format(busted_analyses_dir, ortholog_id)
        os.makedirs(wkd, exist_ok=True)

        # write out the tree to the busted analysis directory in question
        tree_path = '{}/{}_tree.nwk'.format(wkd, ortholog_id)
        with open(tree_path, 'w') as f:
            f.write('{}\n'.format(tree_file))

        # copy the cds fasta alignment from the local directory to this directory

        new_align_path = '{}/{}_qced_cds_aligned.fasta'.format(wkd, ortholog_id)
        shutil.copyfile(orig_align_path, new_align_path)

        # here we now have the alignment and the tree file in the directory we want.

        # now write out the text file that we will use to put in our answeres to the interactive prompts
        answers_script = ['1', '5', '1', new_align_path, tree_path, '1']
        answers_script_path = '{}/{}_answers'.format(wkd, ortholog_id)
        with open(answers_script_path, 'w') as f:
            for line in answers_script:
                f.write('{}\n'.format(line))

        # change dir to the hyphy dir
        output_file_path = '{}/{}_results.out'.format(wkd, ortholog_id)
        os.chdir(HYPHY_dir)
        # HYPHYMP_cmd = './HYPHYMP < {} > {}'.format(answers_script_path, output_file_path)
        # p = subprocess.Popen(HYPHYMP_cmd, shell=True)
        # os.waitpid(p.pid, 0)
        # subprocess.run(['./HYPHYMP', '<', answers_script_path, '>',  output_file_path])
        with open(output_file_path, 'w') as output:
            subprocess.run([HYPHYMP_path], stdout=output, input='\n'.join(answers_script) + '\n', encoding='ascii')

        apples = 'adsf'


    return

def summarise_CODEML_pairwise_analyses():
    ''' Here we will go through each of the block analyses and take out the pairwise dn/ds analyses and put
    them into a single pandas dataframe. The ortholog will the index (which we will evenutally sort) and the
    columns will be each of the pairwise comparisons.'''

    block_analysis_base_dir = '/home/humebc/projects/parky/guidance_analyses'

    list_of_dirs = list()
    for root, dirs, files in os.walk(block_analysis_base_dir):
        list_of_dirs = dirs
        break

    # we will create a index system for the columns so that each pairwise difference will be in a specific column
    # of the df that we are going to create
    # col1 = min_pmin
    # col2 = min_psyg
    # col3 = min_ppsyg
    # col4 = pmin_psyg
    # col5 = pmin_ppsyg
    # col6 = psyg_ppsyg

    comparisons = ['min_pmin', 'min_psyg', 'min_ppsyg', 'pmin_psyg', 'pmin_ppsyg', 'psyg_ppsyg']
    df = pd.DataFrame(columns=comparisons)
    count = 0
    # the list that we will hold the info for a single row in
    dict_of_dnns_values = {}
    for block_dir in list_of_dirs:
        sys.stdout.write('\n\nProcessing block {}\n\n'.format(block_dir))
        # read in the file for this block and populate the df with the dN/dS videos
        with open('{0}/{1}/{1}_guidance_results.out'.format(block_analysis_base_dir, block_dir), 'r') as f:
            out_file = [line.rstrip() for line in f]

        # go line by line through the output file picking up the dnds values and putting them into the df
        for i in range(len(out_file)):
            if 'Data set' in out_file[i]:
                sys.stdout.write('\rProcessing {}'.format(out_file[i]))
            # see if the line in question contains dnds values
            if 'dN/dS=' in out_file[i]:
                # check to see if we have populated a row worth of dnds values
                # if so, populate the df and then start a new row list to collect values in
                if count % 6 == 0 and count != 0:
                    # then we should have a set of dnds values that we can append to the df
                    df.loc[int(ortholog_id)] = [dict_of_dnns_values[pair] for pair in comparisons]
                    dict_of_dnns_values = {}

                # when we are here we are either adding one more value to an already exisiting set
                # or we are starting a brand new set

                # get the dn/ds value and the pair that we are working with
                dn_ds_value = float(out_file[i].split()[7])
                orth_and_pair_info_line = out_file[i-4]
                pair_one = orth_and_pair_info_line.split()[1]
                pair_two = orth_and_pair_info_line.split()[4]
                ortholog_id = pair_one.split('_')[0][1:]
                spp_one = pair_one.split('_')[1][:-1]
                spp_two = pair_two.split('_')[1][:-1]

                if '{}_{}'.format(spp_one, spp_two) in comparisons:
                    dict_of_dnns_values['{}_{}'.format(spp_one, spp_two)] = dn_ds_value
                else:
                    dict_of_dnns_values['{}_{}'.format(spp_two, spp_one)] = dn_ds_value

                count += 1
    # finally add the last set of dn/ds data to the df
    df.loc[int(ortholog_id)] = [dict_of_dnns_values[pair] for pair in comparisons]
    pickle.dump( df, open('{}/CODEML_dnds_pairwise_pandas_df.pickle'.format(block_analysis_base_dir), 'wb'))

    apples = 'asdf'

def generate_summary_stats_for_CODEML_analysis():
    block_analysis_base_dir = '/home/humebc/projects/parky/guidance_analyses'

    df = pickle.load( open('{}/CODEML_dnds_pairwise_pandas_df.pickle'.format(block_analysis_base_dir), 'rb'))
    df.sort_index(inplace=True)
    df.to_csv('/home/humebc/projects/parky/guidance_analyses/CODEML_dN_dS.csv')
    comparisons = ['min_pmin', 'min_psyg', 'min_ppsyg', 'pmin_psyg', 'pmin_ppsyg', 'psyg_ppsyg']

    # # list of orthologs with at least one dnds > 1
    # list_of_orthologs_with_some = []
    # # list of orthologs with only one
    # list_of_orth_with_one = []
    # # list of orthologs with only two
    # list_of_orth_with_two = []
    # # list of orthologs with only three
    # list_of_orth_with_three = []
    # # list of orthologs with only four
    # list_of_orth_with_four = []
    count_dict = defaultdict(list)

    # # list of min_pmin
    # list_of_min_pmin = []
    # # list_of_min_psyg
    # list_of_min_psyg = []
    # # list_of_min_ppsyg
    # list_of_min_ppsyg = []
    # # list of pmin_psyg
    # list_of_pmin_psyg = []
    # # list of pmin_ppsyg
    # list_of_pmin_ppsyg = []
    # #list of psyg_ppsyg
    # list_of_psyg_ppsyg = []
    pairwise_dict = defaultdict(list)

    total_count = 0
    for index in df.index.values.tolist():
        sys.stdout.write('\rProcessing index: {}'.format(index))
        count = 0
        for comp in comparisons:
            if df.loc[index, comp] > 1 and df.loc[index, comp] != float(99):
                count +=1
                total_count += 1
                pairwise_dict[comp].append(index)
        count_dict[count].append(index)

    f = open('{}/CODEML_dN_dS_results.txt'.format(block_analysis_base_dir), 'w')

    print('PAIRWISE COUNTS:')
    f.write('PAIRWISE COUNTS:\n')
    for comp in comparisons:
        print('{}: {}'.format(comp, len(pairwise_dict[comp])))
        f.write('{}: {}\n'.format(comp, len(pairwise_dict[comp])))
    print('\n\nCOUNTS')
    f.write('\n\nCOUNTS\n')
    count = 0
    for i in [1,2,3,4, 5, 6]:
        print('Orthologs with {} ORF dN/dS > 1: {} '.format(i, len(count_dict[i])))
        f.write('Orthologs with {} ORF dN/dS > 1: {} \n'.format(i, len(count_dict[i])))
        count += len(count_dict[i])
    print('\n\nNumber of orthologs with at least one >1 dN/dS scoring ORF: {}'.format(count))
    f.write('\n\nNumber of orthologs with at least one >1 dN/dS scoring ORF: {}\n'.format(count))
    f.close()

    apples = 'asdf'

def summarise_BUSTED_analyses():
    busted_base_dir = '/home/humebc/projects/parky/busted_analyses'

    list_of_dirs = list()
    for root, dirs, files in os.walk(busted_base_dir):
        list_of_dirs = dirs
        break

    manager = Manager()
    managed_list = manager.list()

    input_queue = Queue()

    for dir in list_of_dirs:
        input_queue.put(dir)

    num_proc = 40

    for i in range(num_proc):
        input_queue.put('STOP')

    list_of_processes = []
    for N in range(num_proc):
        p = Process(target=collect_busted_MP_worker, args=(input_queue,managed_list))
        list_of_processes.append(p)
        p.start()

    for p in list_of_processes:
        p.join()




    df = pd.DataFrame([tup[1] for tup in list(managed_list)], columns=['p_value_div_selct'], index=[tup[0] for tup in list(managed_list)])
    # for ortholog_id, p_value_float in list(managed_list):
    #     df.loc[int(ortholog_id)] = p_value_float
    df.sort_index(inplace=True)
    pickle.dump(df, open('{}/BUSTED_p_value_results.pickle'.format(busted_base_dir), 'wb'))
    df.to_csv('{}/BUSTED_p_value_results.csv'.format(busted_base_dir))
    apples = 'asdf'


def collect_busted_MP_worker(input_queue, managed_list):
    busted_base_dir = '/home/humebc/projects/parky/busted_analyses'
    for directory in iter(input_queue.get, 'STOP'):
        sys.stdout.write('\rCollecting ortholog for {}'.format(directory))
        with open('{0}/{1}/{1}_results.out'.format(busted_base_dir, directory)) as f:
            output_file = [line.rstrip() for line in f]

        for i in range(len(output_file)):
            if output_file[i] == 'MATRIX':
                ortholog_id = int(output_file[i + 1].split('_')[0].lstrip()[1:])
            if 'Likelihood ratio test for episodic' in output_file[i]:
                p_value_float = float(output_file[i].split()[-1][:-3])

        managed_list.append((ortholog_id, p_value_float))

def summarise_busted_results_stats():
    busted_base_dir = '/home/humebc/projects/parky/busted_analyses'
    df = pickle.load( open('{}/BUSTED_p_value_results.pickle'.format(busted_base_dir), 'rb'))

    sig_at_dict = defaultdict(list)

    for index in df.index.values.tolist():
        p_val = df.loc[index, 'p_value_div_selct']
        if p_val < 0.05:
            sig_at_dict[0.05].append(index)
        if p_val < 0.01:
            sig_at_dict[0.01].append(index)
        if p_val < 0.001:
            sig_at_dict[0.001].append(index)

    f = open('{}/summary_stats_for_busted_p_vals.txt'.format(busted_base_dir), 'w')

    print('Number of orthologs found to have significant diversifying selection at:\n\n')
    print('p < 0.05 = {}'.format(len(sig_at_dict[0.05])))
    print('p < 0.01 = {}'.format(len(sig_at_dict[0.01])))
    print('p < 0.001 = {}'.format(len(sig_at_dict[0.001])))
    f.write('Number of orthologs found to have significant diversifying selection at:\n\n')
    f.write('p < 0.05 = {}\n'.format(len(sig_at_dict[0.05])))
    f.write('p < 0.01 = {}\n'.format(len(sig_at_dict[0.01])))
    f.write('p < 0.001 = {}\n'.format(len(sig_at_dict[0.001])))
    f.close()

generate_summary_stats_for_CODEML_analysis()


