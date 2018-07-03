import pandas as pd
import os
from multiprocessing import Queue, Process, Manager, current_process
from plumbum import local


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

def convert_interleaved_to_sequencial_fasta(fasta_in):

    list_seq_names = []
    list_seq_sequences = []
    num_seqs = 0

    fasta_cropped = []
    # Get rid of the first line and get rid of the blank lines
    for line in fasta_in[1:]:
        if line != '':
            fasta_cropped.append(line)

    for line in fasta_cropped:
        if line.startswith('>'):
            num_seqs += 1

    for i in range(len(fasta_cropped)):
            if i < num_seqs:
                # Then we are on one of the inital lines
                list_seq_names.append(fasta_cropped[i].split()[0])
                list_seq_sequences.append(''.join(fasta_cropped[i].split()[1:]))
            else:
                index = i%num_seqs
                list_seq_sequences[index] += ''.join(fasta_cropped[i].split()[1:])

    out_fasta = []
    for name, seq in zip(list_seq_names, list_seq_sequences):
        out_fasta.extend(['>{}'.format(name), seq])

    return out_fasta

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


def main():
    # the amino acid sequences
    aa_seq_array = pd.read_csv('/home/humebc/projects/parky/aa_tree_creation/aa_seq.csv', sep=',', lineterminator='\n', index_col=0, header=0)

    # the gene IDs
    gene_id_array = pd.read_csv('/home/humebc/projects/parky/aa_tree_creation/gene_id.csv', sep=',', lineterminator='\n', index_col=0, header=0)

    # the row and column headings should be the same for both arrays so we can iterate both at the same time
    # for each row
    # each row represents an ortholog. These will be our units to process. We may want to multiprocess this.
    # to make this multiprocess friendly, lets through tuples of gene_id_name and aa_sequence into a list and then
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
        p = Process(target=create_local_alignment, args=(ortholog_input_queue, output_dir, spp_list))
        allProcesses.append(p)
        p.start()

    for p in allProcesses:
        p.join()

def create_local_alignment(input, output_dir, spp_list):
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


main()




























