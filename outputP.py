"""This file contains the functions for Diamond output processing"""
__author__ = "6947325: Johannes Zieres"
__credits__ = ""
__email__ = "johannes.zieres@gmail.com"

import os
import glob
import copy
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import pandas as pd
import seaborn as sns


def exclude_putative_transition_frameshift(query_alignment, subject_alignment, frameshift_position):

    """
    Checks if the detected frameshifts correlates with a putative intron an could be thus neglected.
    :param query_alignment: aligned query subsequence
    :param subject_alignment: aligned subject subsequence
    :param frameshift_position: position where a frameshift was detected
    :return: returns True if the frameshift correlates with an putative intron, and False if not
    """

    # checks if the detected frameshift is a result of a putative Intron
    # if case for normal, non edge positions
    if not (((frameshift_position + 18) > len(query_alignment)) or ((frameshift_position - 17) < 0)):
        # later returned, True if the frameshifts correlates with a putative intron and could be neglected
        wrong_intron_left = True
        wrong_intron_right = True

        # used range starts at i-2 to be fault tolerant in intron-exon transition positions
        for y in range(frameshift_position - 17, frameshift_position - 2):  # upstream region
            if subject_alignment[y] != "-":
                wrong_intron_left = False
                break

        for y in range(frameshift_position + 3, frameshift_position + 18):  # downstream region
            if subject_alignment[y] != "-":
                wrong_intron_right = False
                break

        wrong_intron = wrong_intron_left or wrong_intron_right

    # case for frameshifts at the edge of the alignment
    else:
        # stays True if the frameshift might be the result of a putative intron
        wrong_intron = True
        # case for a downstream edge position
        if (frameshift_position + 18) > len(query_alignment):
            for y in range(frameshift_position - 17, frameshift_position - 2):
                if subject_alignment[y] != "-":
                    wrong_intron = False
                    break

        # case for a upstream edge position
        else:
            for y in range(frameshift_position + 3, frameshift_position + 18):
                if subject_alignment[y] != "-":
                    wrong_intron = False
                    break

    return wrong_intron


def read_in_diamond_output(output_dir):
    """
    The function reads in the Diamond blastx results in a given directory. The results for each query are combined and
    stored in the all_diamond_results list, which is also returned. Stores the processed results of each query in the
    form: [scaffold, low_cov_start_pos_in_scaffold, low_cov_end_pos_in_scaffold, query_start_pos_in_scaffold,
    query_end_pos_in_scaffold, Sublist with the information of each query hit: [[protein_hit, e_value, bit_score,
    similarity_percentage, frameshift_list, query_alignment_start_pos, query_alignment_end_pos],
    ...]], which are sublists of the all_diamond_results list. The frameshift list of each protein hit contains all
    detected frameshifts with their positions in the query (I = insertion, D = deletion) e.g. (20, I).
    But not all frameshifts are considered, putative detections at intron/exon transitions are excluded (see
    exclude_putative_transition_frameshift() function)
    PS: frameshift positions are the python list positions.
    :param output_dir: path to the output directory
    :return: all_diamond_results
    """
    # Returns every file in the directory with .out at the end
    output_file_list = glob.glob(output_dir + "/temp_out_*.txt")

    # Sorting the file path ascending by their file number
    dir_path_len = len(output_dir) + 9  # +9 because temp_out_ consists of 9 chars
    output_file_list = sorted(output_file_list, key=lambda current_path: int(current_path[dir_path_len:-4]))

    all_diamond_results = []  # Contains for each query a sublist with all diamond hits

    # Reading out all output files
    for x in output_file_list:
        current_output_file = open(x, "r")

        # Initialising the parameters
        query_hit_list = []    # Contains all hits of the current query

        # Position information of the current region
        # [scaffold, low_cov_start_pos_in_scaffold, low_cov_end_pos_in_scaffold, query_start_pos_in_scaffold,
        # query_end_pos_in_scaffold]
        current_region = ["initscaffold", 0, 0, 0, 0]
        frameshift_list = []  # Contains the positions of detected frameshifts
        query_alignment = ""
        # The following two parameters are important for the localising of the frameshift in the scaffold, and the later
        # Overlapping heuristic for decide which frameshifts are considered in each query
        query_alignment_start_pos = 0   # Saves the start position of the alignment in the query, initialising it with 0
                                        # Is no problem because Diamond starts with counting at 1

        query_alignment_end_pos = 0     # for saving the current end position of each alignment
        subject_alignment = ""
        # The following three store the respective values
        protein_hit = ""
        e_value = ""
        bit_score = ""
        similarity_percentage = ""

        # Reads in the current output file
        for line in current_output_file:

            # Case for a new Query, if for the previous regions hits were found the hits were appended to the outputlist
            if line.startswith("Query="):

                # Case if the query_hit_list, which holds the hits of each query, is not empty,
                # which means the previous region had an Diamond hit
                if query_hit_list:

                    # For every query_hit_list per query, the first initialising element must be deleted and the last
                    # hit must be afterwards be appended

                    # Appends the last remaining protein hit of the query, because each protein hit of the query is
                    # appended to the list by its following hit,which means that the last one needs to be appended later
                    query_hit_list.append([protein_hit, e_value, bit_score, similarity_percentage, query_alignment,
                                           subject_alignment, query_alignment_start_pos, query_alignment_end_pos])

                    # Reset the parameter
                    query_alignment_start_pos = 0

                    # Removes the initialising element of the query, because the first element appends the last protein
                    # hit of the previous query, which is a result of the used approach
                    query_hit_list.pop(0)

                    # Stores the processed results of each query in the form:
                    # [scaffold, low_cov_start_pos_in_scaffold, low_cov_end_pos_in_scaffold,
                    # query_start_pos_in_scaffold, query_end_pos_in_scaffold,
                    # Sublist with the information of each query hit:
                    # [[protein_hit, e_value, bit_score, similarity_percentage, frameshift_list,
                    # query_alignment_start_pos, query_alignment_end_pos], ...]]
                    processed_query_hit_list = [current_region[0], current_region[1], current_region[2],
                                                current_region[3], current_region[4], []]

                    # Processes and appends each protein hit of the query to the processed_query_hit_list and if all
                    # hits are processed, appends those list to the all_diamond_results list, which is the output of
                    # theses function
                    for current_protein_hit in query_hit_list:

                        # !!!!!IMPORTANT STEP!!!!!!
                        # If the gene of the protein hit lays on the reverse strand, the alignment needs to be inverted
                        # to fit for the following processing
                        if current_protein_hit[6] > current_protein_hit[7]:     # Alignment start pos. > end pos.
                            # Swap the positions to fit the inverted alignment
                            current_protein_hit[6], current_protein_hit[7] = current_protein_hit[7], \
                                                                                 current_protein_hit[6]

                            # Reverse the complete query-protein alignment
                            current_protein_hit[4] = current_protein_hit[4][::-1]   # aligned query part
                            current_protein_hit[5] = current_protein_hit[5][::-1]   # aligned protein part

                        # The placed gaps in the query in the alignment are counted for a correct position calculation
                        query_gap_count = 0

                        # Detect frameshifts in aligned query region
                        for alignment_pos, char in enumerate(current_protein_hit[4]):

                            # Counts how many gaps are inserted into the query before the frameshift position
                            if char == "-":
                                query_gap_count += 1

                            if char == "\\":  # Insertion detected
                                # Function returns True if the frameshift correlates with a putative intron
                                wrong_intron = exclude_putative_transition_frameshift(current_protein_hit[4],
                                                                                      current_protein_hit[5],
                                                                                      alignment_pos)

                                if not wrong_intron:
                                    # Saves the positions as the python list positions of the first nucleotide in the
                                    # triplet with the frameshift
                                    frameshift_list.append(((current_protein_hit[6] - 1) +
                                                            (3 * (alignment_pos - query_gap_count)), "I"))

                            elif char == "/":  # Deletion detected
                                # Function returns True if the frameshift correlates with a putative intron
                                wrong_intron = exclude_putative_transition_frameshift(current_protein_hit[4],
                                                                                      current_protein_hit[5],
                                                                                      alignment_pos)

                                if not wrong_intron:
                                    # Saves the positions as the python list positions of the first nucleotide in the
                                    # triplet with the frameshift
                                    frameshift_list.append(((current_protein_hit[6] - 1) +
                                                            (3 * (alignment_pos - query_gap_count)), "D"))

                        # Appends the data of the current protein hit to the query specific list
                        processed_query_hit_list[5].append([current_protein_hit[0], current_protein_hit[1],
                                                            current_protein_hit[2], current_protein_hit[3],
                                                            frameshift_list,
                                                            current_protein_hit[6], current_protein_hit[7]])

                        # Reset values
                        frameshift_list = []

                    # Reset the the list, the if request is inactive till a new protein hit was found and thus
                    # initialising element is created
                    query_hit_list = []

                    # Appends the processed hit results of the current query to the final output list, if the current
                    # query had at least one diamond hit
                    all_diamond_results.append(processed_query_hit_list)

                # Setting the new region
                # Output example "Query= HiC_scaffold_12#1784003#1784396#11111#22222"
                current_line = line.strip().split()[1].split("#")
                # Format: [scaffold, low_cov_start_pos_in_scaffold, low_cov_end_pos_in_scaffold,
                # query_start_pos_in_scaffold, query_end_pos_in_scaffold]
                current_region = current_line[0:5]

            # Line containing the protein id
            if line.startswith(">"):

                # Appends the previous protein hit to the list with all hits of each query
                # If This is the first hit of te current query, the last protein hit of the previous query is appended,
                # but this "initialising" element is later at the Query= if request removed. In a normal case each hit
                # is appended by the following hit of the same query, and the last hit is also appended at the Query= if
                query_hit_list.append([protein_hit, e_value, bit_score, similarity_percentage, query_alignment,
                                      subject_alignment, query_alignment_start_pos, query_alignment_end_pos])

                # Reset/set the new parameters
                protein_hit = line.split()[0][1:]  # containing the protein hit of the current region
                subject_alignment = ""
                query_alignment = ""
                query_alignment_start_pos = 0

            # Line containing score and e-value
            if line.startswith(" Score"):
                current_line = line.strip().split()
                bit_score = current_line[2]
                e_value = current_line[7]

            # Line similarity percentage
            if line.startswith(" Identities"):
                current_line = line.strip().split()
                similarity_percentage = current_line[3][1:-2]

            # Reading in the aligned sequences
            if line.startswith("Query "):
                # Case that a new alignment starts, and the new start position is registered
                if not query_alignment_start_pos:
                    query_alignment_start_pos = int(line.split()[1])

                query_alignment_end_pos = int(line.split()[3])   # Append the currently last alignment position
                query_alignment += line.split()[2]

            # Subject alignment line case
            if line.startswith("Sbjct "):
                subject_alignment += line.split()[2]

        # Appending the hits of the last query if there were hits(same as upper part)
        if query_hit_list:
            # Appends the last remaining protein hit of the query
            query_hit_list.append([protein_hit, e_value, bit_score, similarity_percentage, query_alignment,
                                   subject_alignment, query_alignment_start_pos, query_alignment_end_pos])

            # Removes the initialising element of the query
            query_hit_list.pop(0)

            # Stores the processed results of each query in the form:
            # [scaffold, query_start_pos_in_scaffold, query_end_pos_in_scaffold, low_cov_start_pos_in_scaffold,
            # low_cov_end_pos_in_scaffold, Sublist with the information of each query hit:
            # [[protein_hit, e_value, bit_score, similarity_percentage, query_alignment,
            #  subject_alignment, query_alignment_start_pos, query_alignment_end_pos], ...]]
            processed_query_hit_list = [current_region[0], current_region[1], current_region[2], current_region[3],
                                        current_region[4], []]

            # Processes and appends each protein hit of the query to the processed_query_hit_list and if all hits are
            # processed, appends those list to the all_diamond_results list, which is the output of theses function
            for current_protein_hit in query_hit_list:

                # !!!!!IMPORTANT STEP!!!!!!
                # If the gene of the protein hit lays on the reverse strand, the alignment needs to be inverted
                # to fit for the following processing
                if current_protein_hit[6] > current_protein_hit[7]:  # Alignment start pos. > end pos.
                    # Swap the positions to fit the inverted alignment
                    current_protein_hit[6], current_protein_hit[7] = current_protein_hit[7], \
                                                                     current_protein_hit[6]

                    # Reverse the complete query-protein alignment
                    current_protein_hit[4] = current_protein_hit[4][::-1]  # aligned query part
                    current_protein_hit[5] = current_protein_hit[5][::-1]  # aligned protein part

                # The placed gaps in the query in the alignment are counted for a correct position calculation
                query_gap_count = 0

                # Detect frameshifts in aligned query region
                for alignment_pos, char in enumerate(current_protein_hit[4]):

                    # Counts how many gaps are inserted into the query before the frameshift position
                    if char == "-":
                        query_gap_count += 1

                    if char == "\\":  # Insertion detected
                        # Function returns True if the frameshift correlates with a putative intron
                        wrong_intron = exclude_putative_transition_frameshift(current_protein_hit[4],
                                                                              current_protein_hit[5],
                                                                              alignment_pos)

                        if not wrong_intron:
                            # Saves the positions as the python list positions of the first nucleotide in the
                            # triplet with the frameshift
                            frameshift_list.append(((current_protein_hit[6] - 1) +
                                                    (3 * (alignment_pos - query_gap_count)), "I"))

                    elif char == "/":  # Deletion detected
                        # Function returns True if the frameshift correlates with a putative intron
                        wrong_intron = exclude_putative_transition_frameshift(current_protein_hit[4],
                                                                              current_protein_hit[5],
                                                                              alignment_pos)

                        if not wrong_intron:
                            # Saves the positions as the python list positions of the first nucleotide in the
                            # triplet with the frameshift
                            frameshift_list.append(((current_protein_hit[6] - 1) +
                                                    (3 * (alignment_pos - query_gap_count)), "D"))

                # Appends the data of the current protein hit to the query specific list
                processed_query_hit_list[5].append([current_protein_hit[0], current_protein_hit[1],
                                                   current_protein_hit[2], current_protein_hit[3], frameshift_list,
                                                   current_protein_hit[6], current_protein_hit[7]])

                # Reset the list for the next protein hit
                frameshift_list = []

            # Appends the processed hit results of the current query to the final output list, if the current
            # query had at least one diamond hit
            all_diamond_results.append(processed_query_hit_list)

        current_output_file.close()

    return all_diamond_results


def check_if_frameshift_is_low_cov(max_detect_dist, low_cov_regions, current_scaffold, framesh_pos_in_scaff):
    """
    This functions checks for each frameshift, that was considered by the overlapping heuristic, if the frameshift
    positions really lays in an original detected low coverage region. Which means before merging low coverage regions
    that are close together and the step where low cov regions were expanded to work as queries for Diamond.
    Because only detected frameshifts in low coverage regions are more likely a result of an assembly error. Detected
    frameshifts in normal coverage regions are more likely the result of the detection of a normal inactivating
    mutation.
    To prevent a hard cutoff between regions where frameshifts are considered and regions where frameshift are excluded,
    the max_detect_dist is introduced, which means when this variable is for example 10... Framshifts that are 10 Bp
    next to a low cov region, are still considered in the healing process.

    :param max_detect_dist:
    :param low_cov_regions:
    :param current_scaffold:
    :param framesh_pos_in_scaff:
    :return:
    """

    # Turns True if the detected frameshifts is in a original low cov region
    low_cov_frameshift = False

    # First search matching scaffold
    for scaffold_data in low_cov_regions:

        # Case that matching scaffold was found
        if scaffold_data[0] == current_scaffold:

            for low_cov_region in scaffold_data[1]:

                considered_start_pos = low_cov_region[0] - max_detect_dist
                considered_end_pos = low_cov_region[1] + max_detect_dist

                # Case for finding a low cov region that overlaps the current frameshift position
                if (considered_start_pos <= framesh_pos_in_scaff) and (considered_end_pos >= framesh_pos_in_scaff):
                    low_cov_frameshift = True
                    break

                # Case where all possible low cov regions are passed, and no overlapping could be found
                elif considered_start_pos > framesh_pos_in_scaff:
                    break

            break

    return low_cov_frameshift


def filter_out_relevant_results(all_diamond_results, max_detect_dist, low_cov_regions):
    """
    This functions works as an additional filtering part, where all irrelevant and not project fitting Diamond results
    are filtered out. This means implementing an overlapping heuristic (for more detailed information about this view
    the project wiki), and a filtering step where frameshifts are filtered out when they are not in the low coverage
    region of the query or right next to it. This could happen when a low cov. region is to short and was expanded to
    the minimum base-pair length to work as a query. In this case frameshifts could also be detected at positions of the
    expanded part of the query and not in the original low cov. part. Or when low cov regions were merged together in
    the merge step, and region between them is also part of the query. Those frameshifts contradict to the program
    approach of healing frameshifts at low cov. regions that correlate with bad polishing. To implement a "soft" cutoff
    between considered and "thrown away" frameshift information, frameshifts that are still close to the original low
    cov. region are also considered. "Close" is defined by the max_detect_dist parameter, so eg. max_detect_dist = 10
    -> frameshifts that are 10Bp away from the low cov. region are still considered.
    The first output list (considered_diamond_hits_list), contains all Diamond hits (with the corresponding
    information), that where considered for the healing by the overlapping heuristic (The list format could be seen
    below in the code comments), REMIND all detected frameshifts of each Diamond hit are stored in this list, even when
    they are not considered in the following healing process!!!
    The second output list (healing_region_list), contains the combined frameshift information of all Diamond hits per
    query, which means the scaffold, start position of the query in the scaffold, end position of the query in the
    scaffold, and a list with all positions in the query, where frameshifts are detected and need to be healed later.
    Only queries with considered frameshifts are appended in those list!!!
    But not all frameshifts, of each by the overlapping heuristic considered protein hit, are considered for the healing
    ,like described before, so there could be Diamond hits where only a part of the detected frameshifts are considered.
    The frameshifts are saved as eg. (20, I) or (12, D), which stands for an insertion 20 base-pairs downstream of the
    query start position, or a deletion 12 base-pair downstream.

    The output file considered_diamond_hits_list contains all relevatn information of each in the healing considered
    query in the format like this:
    [[scaffold, low cov. start pos in scaffold, low cov. end pos. ,query start pos., query end pos.,
                # protein_hit, e_value, bit_score, similarity_percentage, [COMPLETEframeshift_list]], ...]

    :param all_diamond_results: Output of the read_in_diamond_output_function
    :param max_detect_dist: The max_detect_distance defines the distance from
                            a detected frameshift position to the original low cov. region, where a frameshift is still
                            considered and not excluded in the further analysis. Values around 10 might be reasonable.
    :param low_cov_regions: output of the read_in_low_cov_tsv_file() function, to heal only in low cov. regions
    :return: considered_diamond_hits_list, healing_region_list
    """

    # healing_region_list contains the final query information for all queries, and functions as output list
    healing_region_list = []
    # Contains the complete information of all for healing considered Diamond hits by the overlapping heuristic,
    # in a format like this:
    # [[scaffold, low cov. start pos in scaffold, low cov. end pos. ,query start pos., query end pos.,
    # protein_hit, e_value, bit_score, similarity_percentage, [COMPLETEframeshift_list]], ...]
    considered_diamond_hits_list = []

    # Counts how many frameshifts are considered in the healing process
    considered_frameshifts_count = 0

    # Counts how many frameshifts are considered by the overlapping heuristic, includes the frameshifts outside of low
    # cov. regions, excludes those at putative intron transition positions.
    frameshifts_detected = 0

    for query_data in all_diamond_results:

        # Saves the regions of the query that are already covered, by saving the start and end positions in the query of
        # all already considered query-protein alignments
        current_query_region_coverage = []

        # Stores all considered frameshifts of each query
        filtered_query_frameshift_list = []

        for diamond_hit in query_data[5]:

            # Initialise the parameters of the current viewed query-protein alignment
            region_overlapping = False  # Saves if the new region is overlapping with a already included region
            alignment_start_pos = diamond_hit[5]
            alignment_end_pos = diamond_hit[6]

            # Checks for all already added regions, if the new region is overlapping with one of them
            # (Overlapping Heuristic)
            for prev_region in current_query_region_coverage:
                if ((alignment_end_pos >= prev_region[0]) and (alignment_start_pos <= prev_region[0])) or \
                        ((alignment_end_pos >= prev_region[1]) and (alignment_start_pos <= prev_region[1])) or \
                        ((alignment_start_pos >= prev_region[0]) and (alignment_end_pos <= prev_region[1])):
                    region_overlapping = True  # overlapping region found
                    break

            # Case were no previous region is overlapping with the current region
            if not region_overlapping:

                # Append the new region to the coverage list
                current_query_region_coverage.append([alignment_start_pos, alignment_end_pos])

                # Stores the frameshifts of the current diamond hit, that were considered in the following healing
                current_diamond_hit_considered_frameshifts = []
                # Checks for each detected frameshift in the protein-query alignment of the current Diamond hit if the
                # frameshift lays in the low cov. region of the current query
                for single_frameshift in diamond_hit[4]:

                    # Appends the frameshift, when it lays into the low cov. region of the query or right next to it (
                    # distance determined by the max_detect_dist). The previous saved position of each frameshift is
                    # based on the start position of the query, not of the low cov. region
                    # query start position and frameshift position
                    frameshift_pos_in_scaff = int(query_data[3]) + single_frameshift[0]

                    # Checks if in general a frameshift is considered by the overlapping heuristic
                    frameshifts_detected += 1

                    current_scaffold = query_data[0]

                    if check_if_frameshift_is_low_cov(max_detect_dist, low_cov_regions, current_scaffold,
                                                      frameshift_pos_in_scaff):
                        # Appends the checked Frameshift to the considered list
                        current_diamond_hit_considered_frameshifts.append((single_frameshift))
                        considered_frameshifts_count += 1

                # Adds the considered frameshifts of the current Diamond hit to the list with all considered frameshifts
                # of the current query
                filtered_query_frameshift_list += current_diamond_hit_considered_frameshifts

                # Add the considered Diamond hit to the list with all considered protein hits and the corresponding data
                # in a format like this:
                # [[scaffold, low cov. start pos in scaffold, low cov. end pos. ,query start pos., query end pos.,
                # protein_hit, e_value, bit_score, similarity_percentage, [COMPLETEframeshift_list]], ...]

                # If considered for the healing, appends the Diamond hit with the considered frameshifts
                if current_diamond_hit_considered_frameshifts:
                    considered_diamond_hits_list.append([query_data[0], query_data[1], query_data[2], query_data[3],
                                                         query_data[4], diamond_hit[0], diamond_hit[1], diamond_hit[2],
                                                         diamond_hit[3], current_diamond_hit_considered_frameshifts])

        # Append the query data to the healing list, if there was an considered frameshift
        # scaffold, query start pos. in scaff., query end pos. in scaff., frameshift correction list
        if filtered_query_frameshift_list:
            healing_region_list.append([query_data[0], query_data[3], query_data[4], filtered_query_frameshift_list])

    print(considered_frameshifts_count, " Frameshifts are considered in the healing process")

    # regardless if they are in the original low cov or nor, putative intron transition frameshift excluded
    print(frameshifts_detected, " Frameshifts are considered by the overlapping heuristic")

    # Calculate the ratio between healed insertions and deletions
    deletion_count = 0
    insertion_count = 0
    for query in healing_region_list:
        for frameshift_pos in query[3]:
            if frameshift_pos[1] == "D":
                deletion_count += 1
            else:
                insertion_count += 1

    print(insertion_count, "Amount of considered insertions in the healing process")
    print(deletion_count, "Amount of considered deletions in the healing process")

    return considered_diamond_hits_list, healing_region_list


def considered_diamond_hit_length_distribution_plot(considered_diamond_hits_list, output_path,
                                                    original_len_distribution):
    """
    This functions simply calculates the length distribution of all low coverage regions that are used  in the healing
    process. And creates an boxplot out of it, that is saved in the handed over directory.
    The displayed percentage in the bars, is the percentage of the original created queries in that length range, that
    were later considered in the healing process/ healed.
    :param considered_diamond_hits_list: output of one of the filter_out_relevant_results function
    :param output_path: complete path to location where the plot should be saved, INCLUDING the plotname and .png
    :param original_len_distribution: length distribution of the low cov. regions that were used as queries
    :return: length_distribution list and  creates a boxplot in the handed over directory
    """

    # Labels for bars
    tick_labels = ["1-5", "6-25", "25-50", "51-100", "101-250", "250-500", "501-1000", "1001-5000", ">5000"]

    length_distribution = [[5, 0], [25, 0], [50, 0], [100, 0], [250, 0], [500, 0], [1000, 0], [5000, 0],
                           [float("inf"), 0]]

    print(considered_diamond_hits_list[0])

    # List which saves the considered low cov regions, to prevent double counting, when more than one protein hit per
    # query was used in the healing step
    considered_low_cov_regions = []

    # Calculating the length distribution
    for protein_hit in considered_diamond_hits_list:
        low_cov_region_pos = [int(protein_hit[1]), int(protein_hit[2])]
        if low_cov_region_pos not in considered_low_cov_regions:
            considered_low_cov_regions.append(low_cov_region_pos)
            region_length = int(low_cov_region_pos[1]) - int(low_cov_region_pos[0])
            for length_range in length_distribution:
                if region_length <= length_range[0]:
                    length_range[1] += 1
                    break

    # Calculate the percentages of how many of the original low cov regions were then considered in the healing
    # regarding to their length class
    percentage_distribution = []

    for position in range(len(length_distribution)):
        considered_percentage = length_distribution[position][1] / (original_len_distribution[position][1]/100)
        percentage_distribution.append(str(round(considered_percentage, 2)) + "%")

    # X-coordinates of left sides of bars
    elements = [x for x in range(1, len(length_distribution) + 1)]

    # Heights of bars
    element_counts = [y[1] for y in length_distribution]

    # Plotting a bar chart
    plt.barh(elements, element_counts, tick_label=tick_labels, height=0.8, color=["yellow"], edgecolor="black")

    percentagelabel_position = max(element_counts) * 0.035

    # Label the bars
    for index, data in enumerate(element_counts):
        plt.text(y=index + 0.9, x=((data/2) - percentagelabel_position), s=f"{percentage_distribution[index]}",
                 fontdict=dict(fontsize=8))

    # Label the picture
    plt.xlabel("Number of low coverage regions", fontsize=11.5)
    plt.ylabel("Length in bp", fontsize=11.5)
    plt.title("Healed low coverage regions length distribution", fontsize=13, fontweight="bold")

    plt.tight_layout()

    print("figure saved")
    # Function saves the plot
    plt.savefig(output_path, dpi=250)

    return length_distribution


def read_in_toga_lossgene_file(unhealed_file_path, healed_file_path, isoforms_file_path, query_annotation_file_path,
                               putative_false_corrected_tsv_path):
    """
    Small function for evaluating the TOGA loss_summ_data.tsv files. It filters out the genes that are marked in the
    Toga run with the unmodified/unhealed as intact and in the healed assembly as lost or missing.
    For each of those genes all isoforms and all location of each isoform is detected and stored with the relating
    position data in the saved output file.
    :param unhealed_file_path: Toga loss_summ_data.tsv file path of unhealed assembly
    :param healed_file_path: Toga loss_summ_data.tsv file path of healed assembly
    :param isoforms_file_path: Toga output file
    :param query_annotation_file_path: Toga annotation file path
    :param putative_false_corrected_tsv_path: path to the output file that will be created
    :return: results_dataframe
    """

    # Read in the tsv files with different separator
    column_names_1 = ["type", "gene_id", "status"]
    dataframe_unhealed = pd.read_csv(unhealed_file_path, sep='\t', header=None, names=column_names_1,
                                     index_col="gene_id")
    dataframe_healed = pd.read_csv(healed_file_path, sep='\t', header=None, names=column_names_1, index_col="gene_id")

    gene_df_unhealed = dataframe_unhealed[dataframe_unhealed['type'] == 'GENE']

    gene_df_healed = dataframe_healed[dataframe_healed['type'] == 'GENE']

    # Create an combined dataframe with the information of the healed and unhealed file
    combined_df = pd.merge(gene_df_unhealed, gene_df_healed, how="inner",
                           suffixes=("_unhealed", "_healed"), on="gene_id").drop(columns=["type_healed", "type_unhealed"])

    # Search for the gene entries, where the gene was intact in the unhealed assembly and lost, uncertain lost or
    # missing in the healed assembly
    putative_false_corrected_df = combined_df[(combined_df['status_unhealed'] == 'I') &
                                              ((combined_df['status_healed'] == 'UL') |
                                               (combined_df['status_healed'] == 'L') |
                                               (combined_df['status_healed'] == 'M'))]

    # Read in the isoforms tsv file
    gene_isoforms_df = pd.read_csv(isoforms_file_path, sep='\t', index_col="gene_id").rename(
        columns={"trascript_id": "transcript_id"})

    # Print value count before first inner merge (original amount of putative wrong healed genes)
    print(putative_false_corrected_df.value_counts())

    # Not necessary, just for checks
    """# Saves the data for later checks
    putative_false_corrected_df.to_csv("/home/johannes/Desktop/trachinus_draco/TOGA_run_1_output/"
                                       "original_putative_wrong_healed.tsv", sep='\t')"""

    # Includes the informations of the isoform-identifieres for each gene
    putative_false_corrected_df = pd.merge(putative_false_corrected_df, gene_isoforms_df, how="inner", on="gene_id")

    # Read in the complete annotation of the assembly
    column_names_2 = ["scaffold", "2", "3", "start_pos", "end_pos", "6", "7", "8", "9"]
    query_annotation_df = pd.read_csv(query_annotation_file_path, sep='\t', header=None, names=column_names_2).drop(
        columns=["2", "6", "7", "8"])

    # Filter out the entries that are not "transcript" in the 3.rd column
    query_annotation_df = query_annotation_df[query_annotation_df["3"] == "transcript"].drop(columns=["3"])
    # Extract the isoform identifier out of the last column and add those as new column
    query_annotation_df["transcript_id"] = query_annotation_df["9"].str.split("\"").str[1]
    # Removes the projection-location specific identifier number at the end of each isoform identifier
    query_annotation_df["transcript_id"] = query_annotation_df["transcript_id"].map(lambda x: x.rstrip("0123456789"
                                                                                                       ).rstrip("."))
    query_annotation_df = query_annotation_df.drop(columns=["9"]).set_index("transcript_id")

    # Print before second inner merge
    print(putative_false_corrected_df)

    # Inner merge to add the positions of each projection for each gene in the assembly
    putative_false_corrected_df = pd.merge(putative_false_corrected_df.reset_index(), query_annotation_df, how="inner",
                                           on="transcript_id")

    # Print final dataframe
    print(putative_false_corrected_df)

    # Print how many of the original regions are left after merging - some could get lost if there are no entries in the
    # query_annotation file
    print("Genes with transcrips in the annotation file:",
          putative_false_corrected_df.groupby("gene_id")["transcript_id"].nunique().count())

    putative_false_corrected_df.to_csv(putative_false_corrected_tsv_path, sep='\t')

    # Result dataframe with overall results
    results_dataframe = combined_df.value_counts()

    return results_dataframe


def check_overlapping_healing_positions_1(putative_wrong_corrected_file_path,
                                          healing_region_list):
    """

    Note! The best hit of the created queries in the database is not in all cases the hit on which the healing decision
    was based, thus the new version 2 of this function was created (searching for the exact hit)

    This function checks each gene that is marked by TOGA in the unhealed assembly as intact and in healed version
    as lost/missing. This means for all locations of all isoforms of all genes is checked, if they are overlapping with
    a healing position, thus false healing might have happened in this region. For each of those detected overlapping,
    a query is created +- 400 Bp around the healing position. Healing positions that are represented multiple times by
    different queries, are filtered out. The start and end position of each query is afterwards appended to the
    filtered_diamond_query_list, which is afterwards returned for the query file creation.
    :param putative_wrong_corrected_file_path: path to the processed toga output of read_in_toga_lossgene_file()
    :param healing_region_list: output of filter_out_relevant_results()
    :return: filtered_diamond_query_list
    """
    # Read in the file
    input_file = open(putative_wrong_corrected_file_path, "r")
    # Skipps the first line
    next(input_file)

    # Stores the information for the query creation (input of the query_files_creation() function)
    diamond_query_list = []

    # Read in current potential wrong corrected region
    for line in input_file:
        splitted_line = line.split()

        current_region = splitted_line[1]
        isoform_id = splitted_line[4]
        current_scaffold = splitted_line[5]
        start_pos_in_scaff = min(int(splitted_line[6]), int(splitted_line[7]))
        end_pos_in_scaff = max(int(splitted_line[6]), int(splitted_line[7]))

        overlapping_found = False

        # Search if there is a overlapping correction position
        for query in healing_region_list:
            # Case for same scaffold
            if current_scaffold == query[0]:
                # Checks for each healing positions, if its overlapping with the current region
                for healing_pos in query[3]:
                    healing_pos_in_scaff = int(query[1]) + int(healing_pos[0])

                    # Case, where healing in the putative wrong healed region was detected
                    if (healing_pos_in_scaff >= start_pos_in_scaff) and (healing_pos_in_scaff <= end_pos_in_scaff):

                        query_start = max(healing_pos_in_scaff - 400, 0)
                        query_end = healing_pos_in_scaff + 400

                        query_header = current_region + "#" + isoform_id + "#" + current_scaffold + "#" + \
                                       str(start_pos_in_scaff) + "#" + str(end_pos_in_scaff) + "#" + \
                                       str(query_start) + "#" + str(query_end)
                        diamond_query_list.append([current_scaffold, [(query_start, query_end,
                                                                       query_header)]])

                        overlapping_found = True
                        break

            if overlapping_found:
                break

    # Contains the filtered output
    filtered_diamond_query_list = []

    # Removes duplicates, when for a single loci multiple isoforms are predicted, thus multiple queries of the same loki
    # would be created, what makes no sense
    loci_list = []
    for query in diamond_query_list:
        already_added = False
        start_pos = query[1][0][0]
        end_pos = query[1][0][1]
        # Checks if the region of the current query is already represented by another query
        for already_added_loci in loci_list:
            if (already_added_loci[0] == start_pos) and (already_added_loci[1] == end_pos):
                already_added = True
                break

        # Appends the query to the loci list
        if not already_added:
            loci_list.append([start_pos, end_pos])
            filtered_diamond_query_list.append(query)

    for x in filtered_diamond_query_list:
        print(x)

    print(len(filtered_diamond_query_list))

    return filtered_diamond_query_list


def check_overlapping_healing_positions_2(putative_wrong_corrected_file_path,
                                          healing_data_tsv_path, output_tsv_path):
    """
    Improved version, that uses the healing_data.tsv file in which for each healing position the underlying Diamond hit
    is saved.

    This function checks each gene that is marked by TOGA in the unhealed assembly as intact and in healed version
    as lost/missing. This means for all locations of all isoforms of all genes is checked, if they are overlapping with
    a healing position in the .tsv, thus false healing might have happened in this region.

    For each of those detected overlapping, a instance is saved in the output file, double representation of healing
    positions is prevented. The output .tsv could then be used to find the healing related Diamond hits in the output
    file.

    :param putative_wrong_corrected_file_path: path to the processed toga output of read_in_toga_lossgene_file()
    :param healing_data_tsv_path: output of create_detailed_healing_information_file()
    :param output_tsv_path: path to the output file, including the file name
    :return: None
    """

    # Read in the healing_data.tsv file
    healing_tsv = open(healing_data_tsv_path, "r")
    healing_postion_data = []    # Stores the data
    # Read in the lines
    for line in healing_tsv:
        healing_postion_data.append(line.strip().split())
    healing_tsv.close()

    # Read in the file with the putative wrong corrected positions
    input_file = open(putative_wrong_corrected_file_path, "r")
    # Skipps the first line
    next(input_file)

    # Stores the data for the ourput
    output_tsv_data = []

    # Read in current potential wrong corrected region
    for isoform_location_data in input_file:
        split_line = isoform_location_data.split()

        current_region = split_line[1]
        isoform_id = split_line[4]
        current_scaffold = split_line[5]
        start_pos_in_scaff = min(int(split_line[6]), int(split_line[7]))
        end_pos_in_scaff = max(int(split_line[6]), int(split_line[7]))

        # Search if there is a overlapping correction position
        for healing_postion in healing_postion_data:
            # Case for same scaffold
            if current_scaffold == healing_postion[0]:

                healing_pos_in_scaff = int(healing_postion[1])

                # Case, where healing in the putative wrong healed region was detected
                if (healing_pos_in_scaff >= start_pos_in_scaff) and (healing_pos_in_scaff <= end_pos_in_scaff):

                    # Check if the found Diamond hit is already represented
                    protein_hit_already_considered = False

                    for entry in output_tsv_data:
                        print(healing_postion[7])
                        print(entry[12])
                        exit()
                        if healing_postion[7] == entry[12]:
                            protein_hit_already_considered = True
                            break

                    if not protein_hit_already_considered:
                        output_tsv_data.append(split_line[1:5] + healing_postion)

                    # Breaks, thus no more entries with different healing positions in the same query are saved
                    break

    input_file.close()

    # Create the output_file
    output_tsv_file = open(output_tsv_path, "w")

    for data_point in output_tsv_data:
        output_tsv_file.write("\t".join(data_point + ["\n"]))

    output_tsv_file.close()

    return None


def search_relating_diamond_alignments(tsv_file_path, diamond_output_dir,  output_file_path):
    """
    This function searches for given entries in the tsv file the matching Diamond alignments, and save the together in
    the output file for easier checking.
    :param tsv_file_path: path to a output file of check_overlapping_healing_positions_2
    :param diamond_output_dir: path to the dir with the diamond output files
    :param output_file_path: path to the output file, including the file name
    :return: None
    """

    # Read in the .tsv file
    tsv_file = open(tsv_file_path, "r")
    tsv_data = []  # Stores the data
    # Read in the lines
    for line in tsv_file:
        tsv_data.append(line.strip().split())
    tsv_file.close()

    # Searches for the corresponding Diamond hit/alignment
    temp_line_storage = []

    return None


def create_toga_result_plot(results_dataframe):
    """
    Simply creates diagram, that displays the "gene-status flow" between the original and the healed assembly TOGA run
    :return:
    """

    print(results_dataframe)

    sns.set_style("whitegrid")

    status_type_list = ["I", "PI", "UL", "L", "M", "PM", "PG"]

    manual_noob_list_total = [60, 2, -61, -8, 9, -2, 0]

    plot_dataframe = pd.DataFrame(data=manual_noob_list_total, index=status_type_list)
    print(plot_dataframe)

    custom_colours = ["forestgreen", "seagreen", "goldenrod", "darkorange", "indianred", "firebrick", "grey"]

    ax_total = plot_dataframe.plot.bar(color=custom_colours, rot=0, figsize=(8, 5))
    ax_total.set_ylim(-130, 70)
    ax_total.legend(loc='upper right')
    ax_total.set_ylabel("amount of genes")
    ax_total.set_title("Total gene status shift")

    figure_total = ax_total.get_figure()
    figure_total.savefig("/home/johannes/Desktop/Toga_total.png", dpi=300)

    manual_noob_list_loss = [[0, 0, -56, -4, -1, 0, 0],
                              [0, 0, -2, 0, 0, 0, 0],
                              [-111, -4, 0, -25, 0, 0, 0],
                              [-10, 0, -21, 0, -9, 0, 0],
                              [0, 0, 0, -1, 0, 0, 0],
                              [0, 0, 0, -2, 0, 0, 0],
                              [0, 0, 0, 0, 0, 0, 0]]

    plot_dataframe = pd.DataFrame(data=manual_noob_list_loss, columns=status_type_list, index=status_type_list)
    print(plot_dataframe)

    custom_colours = ["forestgreen", "seagreen", "goldenrod", "darkorange", "indianred", "firebrick", "grey"]

    ax_total = plot_dataframe.plot.bar(stacked=True, color=custom_colours, rot=0, figsize=(8, 5))
    ax_total.set_ylim(-150, 0)
    ax_total.legend(loc='lower right')
    ax_total.set_ylabel("negative shift")
    ax_total.set_title("Gene status negative shift")

    figure_total = ax_total.get_figure()
    figure_total.savefig("/home/johannes/Desktop/Toga_loss.png", dpi=400)

    manual_noob_list_gain = [[0, 0, 111, 10, 0, 0, 0],
                             [0, 0, 4, 0, 0, 0, 0],
                             [56, 2, 0, 21, 0, 0, 0],
                             [4, 0, 25, 0, 1, 2, 0],
                             [1, 0, 0, 9, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0]]

    plot_dataframe = pd.DataFrame(data=manual_noob_list_gain, columns=status_type_list, index=status_type_list)
    print(plot_dataframe)

    custom_colours = ["forestgreen", "seagreen", "goldenrod", "darkorange", "indianred", "firebrick", "grey"]

    ax_total = plot_dataframe.plot.bar(stacked=True, color=custom_colours, rot=0, figsize=(8, 5), legend=False)
    ax_total.set_ylim(0, 150)
    ax_total.set_ylabel("positive shift")
    ax_total.set_title("Gene status flow")

    figure_total = ax_total.get_figure()
    figure_total.savefig("/home/johannes/Desktop/Toga_gain.png", dpi=400)

    return None


def main():
    print("Output Processing main executed")


if __name__ == '__main__':
    main()
