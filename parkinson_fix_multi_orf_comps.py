import pandas as pd
from collections import defaultdict
import sys
import itertools
from Bio import pairwise2
from multiprocessing import Queue, Manager
import time
from multiprocessing import Queue, Process, Manager, current_process

def do_the_work():
    '''So now we know that there are a few doozies here that we need to take account of.
    1 - the comps were not unique and had sequence variations. These have been made unique in the
    longest_250 versions of the assemblies and so these are the files we should work with in terms of getting DNA
    2 - the comps are not unique across speceis, i.e. each species has a comp00001
    3 - after the above longest_250 processing we can essentially assume that comps are unique within species

    sooo.... some pseudo code to figure this out
    we should work in a dataframe for this
    read in the list of ortholog gene IDs for each species (the comps that make up the 19000) (this is the dataframe)
    then for each species, identify the gene IDs (comps) for which multiple ORFs were made
    go row by row in the dataframe and see if any of the comp IDs (specific to species) are with multiple ORFs
    These are our rows of interest, we need to work within these rows:

    for each comp for each speceis get a list of the orf aa sequences. turn these into a list of lists
    then use itertools.product to get all of the four orfs that could be aligned.
    then for each of these possible combinations do itertools.combinations(seqs, 2) and do
    pairwise distances and get the average score of the pairwise alignments
    we keep the set of four that got the best possible score
    '''

    # First get the ID of the transcript that is related to this Ortholog for each species
    gene_id_df = pd.read_csv('/home/humebc/projects/parky/gene_id.csv', index_col=0)

    # We want to be able to compare how many of the alignments we fixed and how many were OK due to luck
    # To do this we will need a couple of counters but also the currently chosen ORFs
    aa_seq_df = pd.read_csv('/home/humebc/projects/parky/aa_seq_fixed_again.csv', index_col=0)
    multi_orf_counter = 0
    fixed_counter = 0

    # read in the predicted orfs
    # read in the orf aa files for each of the four species
    aa_orf_file_holder_list = []
    for spp in list(gene_id_df):
        with open('/home/baumgas/projects/done/7_John/assemblies_species/annotation/for-dnds/{}_ORFaa.fa'.format(spp),
                  'r') as f:
            aa_orf_file_holder_list.append(convert_interleaved_to_sequencial_fasta_two([line.rstrip() for line in f]))

    comp_to_orfs_dict_holder_list = [defaultdict(list) for spp in list(gene_id_df)]
    orf_to_aa_dict_holder_list = [dict() for spp in list(gene_id_df)]
    for i in range(len(list(gene_id_df))):
        for j in range(len(aa_orf_file_holder_list[i])):
            if aa_orf_file_holder_list[i][j].startswith('>'):
                comp_ID = aa_orf_file_holder_list[i][j].split()[8].split('_')[0]
                orf_ID = aa_orf_file_holder_list[i][j].split()[0][1:]
                comp_to_orfs_dict_holder_list[i][comp_ID].append(orf_ID)
                orf_to_aa_dict_holder_list[i][orf_ID] = aa_orf_file_holder_list[i][j + 1]
        # print out stats for the default dict to see if we agree with Chris at this point
        multi_orf_comps = sum([1 for k, v in comp_to_orfs_dict_holder_list[i].items() if len(v) > 1])
        print('{} comps with > 1 ORF in {}'.format(multi_orf_comps, list(gene_id_df)[i]))

    # here we have the list of dictionaries from which we can get the comp IDs that have multiple ORFs
    # perhaps we shoud convert them to lists now rather than on the fly
    list_of_multi_orf_comps_per_spp = [[k for k, v in comp_to_orfs_dict_holder_list[i].items() if len(v) > 1] for i in range(len(list(gene_id_df)))]

    # now we go row by row through the df and we ask whether any of the comps are in the respective
    # lists of the list_of_multi_orf_comps_per_spp. When we find a row that has at least one in those lists
    # this is a row we will need to do the computation on to figure out which ORFs should be in the alignment
    mp_list = []
    col_labels = list(gene_id_df)
    section_one_time = 0
    section_two_time = 0
    section_three_time = 0


    for index in gene_id_df.index.values.tolist():

        # within each row check to see if the comp is in its respective list
        row_to_check = False
        for i in range(len(col_labels)):
            if gene_id_df.loc[index, col_labels[i]] in list_of_multi_orf_comps_per_spp[i]:
                row_to_check = True
        # here we have checked through each of the comps in the row
        # if row_to_check == True then we need to work through each of the comps and see which ORF combinations
        # produce the best average pairwise distances


        if row_to_check:
            # print('{}\n{}\n{}\n'.format(section_one_time, section_two_time, section_three_time))
            t0 = time.time()
            print('Checking multi-ORF ortholog {}'.format(index), end='\r')
            multi_orf_counter += 1
            list_of_lists_of_possible_orfs = [comp_to_orfs_dict_holder_list[i][gene_id_df.loc[index, col_labels[i]]] for
                                              i in
                                              range(len(col_labels))]
            t1 = time.time()
            section_one_time += t1 - t0
            t2 = time.time()
            # now run itertools.product on this list of lists to get the tuples which are essentially
            # the different orfs that we would be trying to align
            list_of_lists_of_possible_orfs_as_aa = [[] for spp in col_labels]
            for k in range(len(list_of_lists_of_possible_orfs)):
                for orfID in list_of_lists_of_possible_orfs[k]:
                    list_of_lists_of_possible_orfs_as_aa[k].append(orf_to_aa_dict_holder_list[k][orfID])
            # hardcode it to a list so that we can get the index of each tuple below
            alignment_tuples = [tup for tup in itertools.product(*list_of_lists_of_possible_orfs_as_aa)]
            t3 = time.time()
            section_two_time += t3-t2

            mp_list.append((alignment_tuples, aa_seq_df.loc[index].tolist(), index))




    num_proc = 20

    #Queue that will hold the index of the rows that need to be checked
    input_queue = Queue()

    # populate input_queue
    for tup in mp_list:
        input_queue.put(tup)

    for i in range(num_proc):
        input_queue.put('STOP')

    # Manager for a dict rows_to_be_replaced_dict that will hold the new aa_seqs for the fixed indices
    manager = Manager()
    rows_to_be_replaced_dict = manager.dict()

    list_of_processes = []
    for N in range(num_proc):
        p = Process(target=ORF_screening_worker, args=(input_queue, rows_to_be_replaced_dict))
        list_of_processes.append(p)
        p.start()

    for p in list_of_processes:
        p.join()

    fixed_counter = len(rows_to_be_replaced_dict.keys())
    # now replace the dataframe values
    for index in aa_seq_df.index.values.tolist():
        if index in rows_to_be_replaced_dict.keys():
            # then this is a row that needs replcing with the new values
            aa_seq_df.loc[index] = rows_to_be_replaced_dict[index]


    # at this point it only remains to write the data frame out as csv
    print('{} orthologs were checked due to multi-ORFs\n'
          '{} were fixed\n'
          '{} already contained the optimal choice of ORFs\n'.format
        (
        multi_orf_counter, fixed_counter, multi_orf_counter-fixed_counter
        )
    )
    aa_seq_df.to_csv('/home/humebc/projects/parky/aa_seq_multi_orf_orths_fixed.csv')

    # here we should have a dict with the chosen sequences in and that's really what we were after.
    apples = 'asdf'


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

def ORF_screening_worker(input_queue, rows_to_be_replaced_dict):
    for tup in iter(input_queue.get, 'STOP'):
        alignment_tuples, current_seqs, index = tup
        print('Processing multi-ORF ortholog: {}'.format(index))
        # we should now have sets that are four aa sequences.
        # for each set, generate pairwise comparisons and calculate pairwise distances
        # keep track of each distance and calculate average PW distances for each set_of_alignemtn_seqs
        average_distance_list = []
        for set_of_alignment_seqs in alignment_tuples:
            temp_pairwise_scores_list = []
            for seq_a, seq_b in itertools.combinations(set_of_alignment_seqs, 2):
                # here is a single PW distance calculation
                score = pairwise2.align.globalxx(seq_a, seq_b, score_only=True)
                temp_pairwise_scores_list.append(score)
            # now calculate average
            temp_average = sum(temp_pairwise_scores_list) / len(temp_pairwise_scores_list)
            average_distance_list.append(temp_average)

        # now we have a list of PW distance for each of the sets of sequences (virtual alignments)
        # we want to select the set that has the highest score

        index_of_best_set = average_distance_list.index(max(average_distance_list))
        # now check to see if these match those that were already chosen
        best_set_of_aas = alignment_tuples[index_of_best_set]
        alignments_are_same = True

        for i in range(len(current_seqs)):
            if current_seqs[i] != best_set_of_aas[i]:
                alignments_are_same = False
                break

        if not alignments_are_same:
            # print('Fasta for ortholog: {}'.format(index))
            # for n in range(len(best_set_of_aas)):
            #     print('>{}\n{}'.format(n, alignment_tuples[index_of_best_set][n]))
            # we don't need to know the actual ORFs that have been chosen for each spp
            # only the aa codes so let's just output this result
            rows_to_be_replaced_dict[index] = [aa for aa in alignment_tuples[index_of_best_set]]



do_the_work()