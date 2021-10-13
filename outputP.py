"""This file contains the functions for Diamond output processing"""
__author__ = "6947325: Johannes Zieres"
__credits__ = ""
__email__ = "johannes.zieres@gmail.com"

import os
import glob
import copy
import matplotlib.pyplot as plt


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


def filter_out_relevant_results(all_diamond_results, max_detect_dist):
    """
    This functions works as an additional filtering part, where all irrelevant and not project fitting Diamond results
    are filtered out. This means implementing an overlapping heuristic (for more detailed information about this view
    the project wiki), and a filtering step where frameshifts are filtered out when they are not in the low coverage
    region of the query or right next to it. This could happen when a low cov. region is to short and was expanded to
    the minimum base-pair length to work as a query. In this case frameshifts could also be detected at positions of the
    expanded part of the query and not in the original low cov. part. Those frameshifts contradict to the program
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
    :param all_diamond_results: Output of the read_in_diamond_output_function
    :param max_detect_dist: The max_detect_distance defines the distance from
                            a detected frameshift position to the original low cov. region, where a frameshift is still
                            considered and not excluded in the further analysis. Values around 10 might be reasonable.
    :return: considered_diamond_hits_list, healing_region_list
    """

    # healing_region_list contains the final query information for all queries, and functions as output list
    healing_region_list = []
    # Contains the complete information of all for healing considered Diamond hits by the overlapping heuristic,
    # in a format like this:
    # [[scaffold, low cov. start pos in scaffold, low cov. end pos. ,query start pos., query end pos.,
    # protein_hit, e_value, bit_score, similarity_percentage, [COMPLETEframeshift_list]], ...]
    considered_diamond_hits_list = []

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

                    # Checks the overlapping
                    # (Exlude frameshifts that are not in the low cov. region, or close to them)
                    if ((frameshift_pos_in_scaff >= (int(query_data[1]) - max_detect_dist)) and
                       (frameshift_pos_in_scaff <= (int(query_data[2]) + max_detect_dist))):

                        # Appends the checked Frameshift to the
                        current_diamond_hit_considered_frameshifts.append((single_frameshift))

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

    return considered_diamond_hits_list, healing_region_list


def considered_diamond_hit_length_distribution_plot(considered_diamond_hits_list, output_path):
    """
    This functions simply calculates the length distribution of all low coverage regions that are used  in the healing
    process. And creates an boxplot out of it, that is saved in the handed over directory.
    :param considered_diamond_hits_list: output of one of the filter_out_relevant_results function
    :param output_path: complete path to location where the plot should be saved, INCLUDING the plotname and .png
    :return: length_distribution list and  creates a boxplot in the handed over directory
    """

    # Labels for bars
    tick_labels = ["0-5", "6-25", "25-50", "51-100", "101-250", "250-500", "501-1000", "1001-5000", ">5000"]

    length_distribution = [[5, 0], [25, 0], [50, 0], [100, 0], [250, 0], [500, 0], [1000, 0], [5000, 0],
                           [float("inf"), 0]]

    # Calculating the length distribution
    for protein_hit in considered_diamond_hits_list:
        region_length = int(protein_hit[2]) - int(protein_hit[1])

        for length_range in length_distribution:
            if region_length <= length_range[0]:
                length_range[1] += 1
                break

    # X-coordinates of left sides of bars
    elements = [x for x in range(1, len(length_distribution) + 1)]

    # Heights of bars
    element_counts = [y[1] for y in length_distribution]

    # Plotting a bar chart
    plt.barh(elements, element_counts, tick_label=tick_labels, height=0.8, color=["limegreen"], edgecolor="black")

    # Label the picture
    plt.xlabel("Number of low coverage regions", fontsize=11.5)
    plt.ylabel("Length in bp", fontsize=11.5)
    plt.title("Low coverage regions length distribution", fontsize=13, fontweight="bold")
    plt.xticks(ticks=[20000, 40000, 60000, 80000, 100000, 120000, 140000])
    plt.xlim(xmax=145000)

    plt.tight_layout()

    # Function saves the plot
    plt.savefig(output_path, dpi=250)

    return length_distribution


def main():
    print("Output Processing main executed")


if __name__ == '__main__':
    main()
