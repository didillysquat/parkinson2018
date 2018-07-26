import pandas as pd
from collections import defaultdict
import sys
import itertools
from Bio import pairwise2
from multiprocessing import Queue, Manager
import time
from multiprocessing import Queue, Process, Manager, current_process
import pickle

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
    base_dir = '/home/humebc/projects/parky'
    gene_id_df = pd.read_csv('/home/humebc/projects/parky/gene_id.csv', index_col=0)



    # We want to be able to compare how many of the alignments we fixed and how many were OK due to luck
    # To do this we will need a couple of counters but also the currently chosen ORFs
    aa_seq_df = pd.read_csv('/home/humebc/projects/parky/aa_seq.csv', index_col=0)

    # we will also need to update the cds_seq_df dataframe
    cds_seq_df = pd.read_csv('/home/humebc/projects/parky/cds.csv', index_col=0)


    multi_orf_counter = 0


    # read in the predicted orfs
    # read in the orf aa files for each of the four species
    try:
        comp_to_orfs_dict_holder_list = pickle.load( open('{}/comp_to_orfs_dict_holder_list.pickled'.format(base_dir), 'rb'))
        orf_to_aa_dict_holder_list = pickle.load( open('{}/orf_to_aa_dict_holder_list.pickled'.format(base_dir), 'rb'))
    except:
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

        pickle.dump( comp_to_orfs_dict_holder_list, open('{}/comp_to_orfs_dict_holder_list.pickled'.format(base_dir), 'wb'))
        pickle.dump(orf_to_aa_dict_holder_list, open('{}/orf_to_aa_dict_holder_list.pickled'.format(base_dir), 'wb'))

    # Here we will do the same as above but for the codon sequences.
    # we will create both the com to orf dict and the orf to cds dict.
    # in theory the comp to orf dict should be eactly the same for the cds file as it was for the aa file
    # we should check this as a sanity check. If it is the same then we can ignore this second dict we made

    try:
        comp_to_orfs_dict_holder_list_cds = pickle.load( open('{}/comp_to_orfs_dict_holder_list_cds.pickled'.format(base_dir), 'rb'))
        orf_to_cds_dict_holder_list = pickle.load( open('{}/orf_to_cds_dict_holder_list.pickled'.format(base_dir), 'rb'))
    except:
        cds_orf_file_holder_list = []
        for spp in list(gene_id_df):
            with open('/home/baumgas/projects/done/7_John/assemblies_species/annotation/for-dnds/cds/{}_ORFcds.fa'.format(spp),
                      'r') as f:
                cds_orf_file_holder_list.append(convert_interleaved_to_sequencial_fasta_two([line.rstrip() for line in f]))

        comp_to_orfs_dict_holder_list_cds = [defaultdict(list) for spp in list(gene_id_df)]
        orf_to_cds_dict_holder_list = [dict() for spp in list(gene_id_df)]
        for i in range(len(list(gene_id_df))):
            for j in range(len(cds_orf_file_holder_list[i])):
                if cds_orf_file_holder_list[i][j].startswith('>'):
                    comp_ID = cds_orf_file_holder_list[i][j].split()[8].split('_')[0]
                    orf_ID = cds_orf_file_holder_list[i][j].split()[0][1:]
                    comp_to_orfs_dict_holder_list_cds[i][comp_ID].append(orf_ID)
                    orf_to_cds_dict_holder_list[i][orf_ID] = cds_orf_file_holder_list[i][j + 1]
            # print out stats for the default dict to see if we agree with Chris at this point
            multi_orf_comps = sum([1 for k, v in comp_to_orfs_dict_holder_list_cds[i].items() if len(v) > 1])
            print('{} comps with > 1 ORF in {}'.format(multi_orf_comps, list(gene_id_df)[i]))

        pickle.dump(comp_to_orfs_dict_holder_list_cds,
                    open('{}/comp_to_orfs_dict_holder_list_cds.pickled'.format(base_dir), 'wb'))
        pickle.dump(orf_to_cds_dict_holder_list, open('{}/orf_to_cds_dict_holder_list.pickled'.format(base_dir), 'wb'))
    # SANITY CHECK: PASSED!
    # let's look to make sure that both of the comp_to_orfs dicts are the same from reading in the aa and cds files.
    # for i in range(len(list(gene_id_df))):
    #     comp_dict_aa = comp_to_orfs_dict_holder_list[i]
    #     comp_dict_cds = comp_to_orfs_dict_holder_list_cds[i]
    #     for comp_key in comp_dict_aa.keys():
    #         if set(comp_dict_aa[comp_key]) != set(comp_dict_cds[comp_key]):
    #             sys.exit('The list of ORFs for comp {} do not match between the cds and aa dict'.format(comp_key))


    # here we have the list of dictionaries from which we can get the comp IDs that have multiple ORFs
    # perhaps we shoud convert them to lists now rather than on the fly
    list_of_multi_orf_comps_per_spp = [[k for k, v in comp_to_orfs_dict_holder_list[i].items() if len(v) > 1] for i in range(len(list(gene_id_df)))]

    # now we go row by row through the df and we ask whether any of the comps are in the respective
    # lists of the list_of_multi_orf_comps_per_spp. When we find a row that has at least one in those lists
    # this is a row we will need to do the computation on to figure out which ORFs should be in the alignment
    mp_list = []


    #  N.B. we have a problem here in that some of the CDS sequences are not 3%0.
    # On first inspection it appears that 83 of the orthologs have CDS seqs that are not 3%0.
    # 77 of these happen to be comps that have multi ORFs with 6 not having multi-ORFs
    # I think it is very dangerous to be dealing with these partial ORFs as we don't know if they are partial in the
    # 3' end or the 5' end.
    # As such, I want to get rid of them using the following logic
    # If one of an orthologs cds seqs is not 3%0 and doesn't have multi ORFs then we just drop this ortholog from
    # the study.
    # If above, but there is multi-ORFs then we will first find the best set of ORFs and then check to see if they
    # contain any non 3%0 cds seqs. If they do then we will drop this ortholog too.

    col_labels = list(gene_id_df)
    ortholog_indices_to_drop = []
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
            sys.stdout.write('\rChecking multi-ORF ortholog {} for MP'.format(index))
            multi_orf_counter += 1
            list_of_lists_of_possible_orfs = [comp_to_orfs_dict_holder_list[i][gene_id_df.loc[index, col_labels[i]]] for
                                              i in
                                              range(len(col_labels))]


            # convert the orf IDs associated to each of the possible ORFs to both aa_seqs and cds_seqs
            list_of_lists_of_possible_orfs_as_aa = [[] for spp in col_labels]
            list_of_lists_of_possible_orfs_as_cds = [[] for spp in col_labels]
            for k in range(len(list_of_lists_of_possible_orfs)):
                for orfID in list_of_lists_of_possible_orfs[k]:
                    list_of_lists_of_possible_orfs_as_aa[k].append(orf_to_aa_dict_holder_list[k][orfID])
                    list_of_lists_of_possible_orfs_as_cds[k].append(orf_to_cds_dict_holder_list[k][orfID])


            # hardcode it to a list so that we can get the index of each tuple below

            # to ensure maintainance of the direct link between the aa seq and the cds_seq it is maybe easiest
            # if instead of the unit of the tup being an aa sequence, it is a tuple itself which is (aa_seq, cds_seq)
            # this way we can know which of the cds_seqs were chosen based on the AAs/ORFIDs.
            # Create this list of list of tuples using the list_of_lists_of_possible_orfs and ..._as_aa lists
            list_of_lists_of_aa_to_cds_tups = [[] for spp in col_labels]
            for n in range(len(list_of_lists_of_possible_orfs)):
                for m in range(len(list_of_lists_of_possible_orfs[n])):
                    aa_seq = list_of_lists_of_possible_orfs_as_aa[n][m]
                    cds_seq = list_of_lists_of_possible_orfs_as_cds[n][m]
                    list_of_lists_of_aa_to_cds_tups[n].append((aa_seq, cds_seq))

            # With the new change we should now be passing through a set of tuples that contain the aa_seq and the
            # ORF id to the MP
            alignment_tuples = [tup for tup in itertools.product(*list_of_lists_of_aa_to_cds_tups)]


            mp_list.append((alignment_tuples, aa_seq_df.loc[index].tolist(), cds_seq_df.loc[index].tolist(), index))
        else:
            # If the row does not need checking for multi-ORFs then we should still check to make sure that
            # the cds seqs for this ortholog are modulus 3 compliant else we should get rid of the ORF
            drop = False
            for spp in col_labels:
                if len(cds_seq_df.loc[index, spp])%3 != 0:
                    drop = True
                    break
            if drop:
                ortholog_indices_to_drop.append(index)


    #here we need to make some changes. The dN/dS analysis that we are going to perform will be based on
    # codon-based DNA sequences rather than amino acid seuences. As such we will need to have some way of
    # correcting the cds dataframe as well as the AA dataframe.
    num_proc = 20

    #Queue that will hold the index of the rows that need to be checked
    input_queue = Queue()

    # populate input_queue
    for tup in mp_list:
        input_queue.put(tup)

    for i in range(num_proc):
        input_queue.put('STOP')

    # Manager for a dict rows_to_be_replaced_dict that will hold the new (aa_seqs, ORF_id) for the fixed indices
    manager = Manager()
    rows_to_be_replaced_dict = manager.dict()
    list_of_indices_to_drop = manager.list()

    list_of_processes = []
    for N in range(num_proc):
        p = Process(target=ORF_screening_worker, args=(input_queue, rows_to_be_replaced_dict, list_of_indices_to_drop))
        list_of_processes.append(p)
        p.start()

    for p in list_of_processes:
        p.join()

    # add the multi proc list of indices to drop to the non-multi proc list
    ortholog_indices_to_drop.extend(list(list_of_indices_to_drop))

    # print out the orthologs that were non %3 compliant
    with open('non_mod_3_compliant_orthologs.txt', 'w') as f:
        for line in ortholog_indices_to_drop:
            f.write('{}\n'.format(line))

    # print out the orthologs that needed fixing
    with open('multi_orf_orthologs_that_needed_fixing.txt', 'w') as f:
        for orth_id in rows_to_be_replaced_dict.keys():
            f.write('{}\n'.format(orth_id))

    # now drop the indices from the dfs (all three)
    aa_seq_df = aa_seq_df.drop(ortholog_indices_to_drop, axis=0)
    cds_seq_df = cds_seq_df.drop(ortholog_indices_to_drop, axis=0)
    gene_id_df = gene_id_df.drop(ortholog_indices_to_drop, axis=0)

    fixed_counter = len(rows_to_be_replaced_dict.keys())
    # now replace the dataframe values
    # we can use a single parse on both the aa and cds
    for index in aa_seq_df.index.values.tolist():
        if index in rows_to_be_replaced_dict.keys():
            # then this is a row that needs replcing with the new values
            aa_seq_df.loc[index] = [tup[0] for tup in rows_to_be_replaced_dict[index]]
            cds_seq_df.loc[index] = [tup[1] for tup in rows_to_be_replaced_dict[index]]


    # at this point it only remains to write the dataframes out as csv
    print('{} orthologs were checked due to multi-ORFs\n'
          '{} were fixed\n'
          '{} already contained the optimal choice of ORFs\n'.format
        (
        multi_orf_counter, fixed_counter, multi_orf_counter-fixed_counter
        )
    )
    aa_seq_df.to_csv('/home/humebc/projects/parky/aa_seq_multi_orf_orths__partial_cdsfixed.csv')
    cds_seq_df.to_csv('/home/humebc/projects/parky/cds_seq_multi_orf_orths__partial_cdsfixed.csv')
    gene_id_df.to_csv('/home/humebc/projects/parky/gene_id_multi_orf_orths_partial_cds_fixed.csv')

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

def ORF_screening_worker(input_queue, rows_to_be_replaced_dict, indices_to_drop_list):
    #also verify that the cds sequences chosen are the same or different
    # to do this we will have to pass in the current cds seqs as well.

    for tup in iter(input_queue.get, 'STOP'):
        alignment_tuples, current_seqs_aa, current_seqs_cds, index = tup
        sys.stdout.write('\rAssessing multi-ORF ortholog: {} for best alignments'.format(index))

        # we should now have sets that are four aa sequences each paired with four cds seqs.
        # for each set, generate pairwise comparisons and calculate pairwise distances based on the aa_seqs
        # keep track of each distance and calculate average PW distances for each set_of_alignemtn_seqs
        # also associate which set of ORF_ids produced each average PW distance.
        # This list will keep track of this by storing tups that are (average PW distance)
        average_distance_list = []
        # it is important to note, and maintain the fact that the aa seqs and ORF_ids are in the species order
        # inside the alignment_tupes, i.e. first seq is always min
        for set_of_alignment_seqs in alignment_tuples:
            temp_pairwise_scores_list = []
            # here we are comparing two tuples that are (aa_seq, ORF_id)

            for seq_a, seq_b in itertools.combinations(set_of_alignment_seqs, 2):

                # here is a single PW distance calculation
                score = pairwise2.align.globalxx(seq_a[0], seq_b[0], score_only=True)
                temp_pairwise_scores_list.append(score)
            # now calculate average
            temp_average = sum(temp_pairwise_scores_list) / len(temp_pairwise_scores_list)
            average_distance_list.append(temp_average)

        # now we have a list of PW distance for each of the sets of sequences (virtual alignments)
        # we want to select the set that has the highest score

        index_of_best_set = average_distance_list.index(max(average_distance_list))
        # now check to see if these match those that were already chosen
        best_set_of_aa = [tup[0] for tup in alignment_tuples[index_of_best_set]]
        best_set_of_cds = [tup[1] for tup in alignment_tuples[index_of_best_set]]
        alignments_are_same = True

        for i in range(len(current_seqs_aa)):
            if current_seqs_aa[i] != best_set_of_aa[i]:
                alignments_are_same = False
                break

        if alignments_are_same:
            for i in range(len(current_seqs_cds)):
                if current_seqs_cds[i] != best_set_of_cds[i]:
                    alignments_are_same = False
                    break
        for i in range(len(current_seqs_cds)):
            if len(best_set_of_cds[i])%3 != 0:
                # if any of the new cds seqs break the 3%0 rule then we don't want to drop the ortholog
                indices_to_drop_list.append(index)
                # set alignments to true so that this ortholog will not be modified in the df s
                alignments_are_same = True
                break

        if not alignments_are_same:

            # we don't need to know the actual ORFs that have been chosen for each spp
            # only the aa codes so let's just output this result
            rows_to_be_replaced_dict[index] = [aa_orf_id_tup for aa_orf_id_tup in alignment_tuples[index_of_best_set]]



do_the_work()