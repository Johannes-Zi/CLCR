"""This program is for the database comparison of the detected regions"""
__author__ = "6947325: Johannes Zieres"
__credits__ = ""
__email__ = "johannes.zieres@gmail.com"

import os
import glob
import copy


def database_comparison(list_with_scaffold_specific_low_cov_reg_lists, fna_file_path, min_length, query_dir,
                        seq_per_fasta):

    """
    This function creates the .fasta query files for a diamond slurm job, out of a input list with the low coverage
    regions of each and a path to the matching fna file to read out the sequences of the low coverage regions.
    The saved positions of the queries in each scaffold are the normal positions, not python list positions(-1) !!!
    :param list_with_scaffold_specific_low_cov_reg_lists: list consisting of scaffold specific sublists with tuples
                                                          containing the start and he end of a region for database
                                                          comparison (output of detect_regions, merge_regions or
                                                          combine_regions_multiple_scaffolds).
    :param fna_file_path: path of the fna file for receiving the base sequence
    :param min_length: expanding the detected regions to 500 if they are smaller, for a better working on diamond
    :param query_dir: directory where the query files are stored in
    :param seq_per_fasta: max amount of queries stored per fasta, The fisrt and the last file could contain in special
                          cases more queries
    :return no return
    """

    # Read in the .fna file for fast access
    input_fna_file = open(fna_file_path)
    scaffold_list = []  # Filled with the sequences of the fna file
    temp_scaffold = []  # Contains the current scaffold
    current_header = ""  # Containing the header of the current scaffold

    # Filling the .fna region list/ reading in the original assembly
    for line in input_fna_file:

        # Sequence line case
        if (line[0] != ">") and (line[0] != ";"):
            for x in line.strip():  # filling up a new sequence
                temp_scaffold.append(x)

        # Sequence header case
        elif line[0] == ">":
            scaffold_list.append([current_header, temp_scaffold])  # appending the previous region
            current_header = line.strip().split()[0][1:]  # reading in the new header
            temp_scaffold = []  # resetting the list/ ready for new region

    # Append last remaining region
    scaffold_list.append([current_header, temp_scaffold])
    input_fna_file.close()
    # Delete the first empty initialising element
    scaffold_list.pop(0)

    # Initialise some parameters
    fasta_count = 0     # for individual naming at creating of the .fasta files
    count_added_queries = 0     # Counts how many queries are created overall
    already_added_queries = 0   # Counts how many queries are already added to the fasta file

    # Creating the first fasta file
    new_fasta = query_dir + "/" + "temp_in_" + str(fasta_count) + ".fasta"
    current_fasta = open(new_fasta, "w")

    # Creating the query files for all scaffolds, searching a the matching scaffold sequence for each
    # current_low_cov_list with the low cov. regions of each scaffold
    for scaffold in list_with_scaffold_specific_low_cov_reg_lists:
        scaffold_name = scaffold[0]
        current_low_cov_list = scaffold[1]
        current_tuple = 0  # saves the position of the current region tuple in current_low_cov_list

        for current_scaffold in scaffold_list:
            if scaffold_name == current_scaffold[0]:

                # Possible first cases where the region cant be expanded 250 to the left,
                # because the (startposition(of the region) - 250) < 0
                while (current_low_cov_list[current_tuple][0] - int(min_length / 2)) <= 0:

                    # Calculating the expand length
                    region_length = current_low_cov_list[current_tuple][1] - current_low_cov_list[current_tuple][0]
                    if region_length < min_length:
                        expand_left = expand_right = int((min_length - region_length) / 2)
                    else:
                        expand_left = expand_right = 0

                    # Case if the starting point is moved in the negative by expand_left
                    if (current_low_cov_list[current_tuple][0] - expand_left) < 0:
                        region_start = 0
                        region_end = current_low_cov_list[current_tuple][1] + expand_right + \
                                     abs(current_low_cov_list[current_tuple][0] - expand_left)
                    # Normal case
                    else:
                        region_start = current_low_cov_list[current_tuple][0] - expand_left
                        region_end = current_low_cov_list[current_tuple][1] + expand_right

                    # Writing the region in the fasta file
                    current_fasta.write((">" + scaffold_name + "#" + str(current_low_cov_list[current_tuple][0]) +
                                         "#" + str(current_low_cov_list[current_tuple][1]) + "\n"))
                    current_fasta.write(("".join(current_scaffold[1][region_start:region_end]) + "\n"))
                    current_tuple += 1
                    count_added_queries += 1
                    already_added_queries += 1

                    # Break condition
                    if current_tuple >= len(current_low_cov_list):
                        break

                break_value = False
                # Each while run fills up a new fasta file
                while True:

                    # Break condition
                    if (current_tuple >= len(current_low_cov_list)) or \
                            ((current_low_cov_list[current_tuple][1]) >
                             (len(current_scaffold[1]) - int(min_length / 2))):
                        break_value = True
                        break

                    # Writes the header
                    current_fasta.write((">" + scaffold_name + "#" + str(current_low_cov_list[current_tuple][0]) +
                                         "#" + str(current_low_cov_list[current_tuple][1]) + "\n"))

                    # Calculation the length a region must be expanded
                    region_length = current_low_cov_list[current_tuple][1] - current_low_cov_list[current_tuple][0]

                    if region_length < min_length:
                        expand_len = int((min_length - region_length) / 2)
                    else:
                        expand_len = 0

                    # Write the query sequence in the file
                    current_fasta.write("".join(current_scaffold[1]
                                                [(current_low_cov_list[current_tuple][0] - expand_len):
                                                 (current_low_cov_list[current_tuple][1] + expand_len)])
                                        + "\n")

                    # Update parameters
                    current_tuple += 1
                    count_added_queries += 1
                    already_added_queries += 1

                    # Break condition
                    if break_value:
                        break

                    # Creating .fasta file, if the current file is full (defined by seq_per_fasta)
                    if already_added_queries >= seq_per_fasta:
                        fasta_count += 1
                        current_fasta.close()
                        already_added_queries = 0
                        new_fasta = query_dir + "/" + "temp_in_" + str(fasta_count) + ".fasta"
                        current_fasta = open(new_fasta, "w")

                # Possible last cases where the region cant be expanded 250 to the right,
                # because (the endposition(of the region) + 250) > sequence length
                while current_tuple < len(current_low_cov_list):
                    # Calculating the expand length
                    region_length = current_low_cov_list[current_tuple][1] - current_low_cov_list[current_tuple][0]
                    if region_length < min_length:
                        expand_left = expand_right = int((min_length - region_length) / 2)
                    else:
                        expand_left = expand_right = 0

                    # Case if the endpoint is moved behind the end of the current_scaffold by expand_right
                    if (current_low_cov_list[current_tuple][1] + expand_right) > len(current_scaffold[1]):
                        region_start = current_low_cov_list[current_tuple][0] - expand_left - \
                                       ((current_low_cov_list[current_tuple][1] + expand_right) -
                                        len(current_scaffold[1]))
                        region_end = len(current_scaffold[1])
                    # Normal case
                    else:
                        region_start = current_low_cov_list[current_tuple][0] - expand_left
                        region_end = current_low_cov_list[current_tuple][1] + expand_right

                    # Writing the region in the fasta file
                    current_fasta.write((">" + scaffold_name + "#" + str(current_low_cov_list[current_tuple][0]) +
                                         "#" + str(current_low_cov_list[current_tuple][1]) + "\n"))
                    current_fasta.write(("".join(current_scaffold[1][region_start:region_end]) + "\n"))
                    current_tuple += 1
                    count_added_queries += 1

                break       # Break because there is only one matching scaffold sequence for each scaffold

    # Close the last .fasta file
    current_fasta.close()
    fasta_count += 1

    # User information
    print("Amount of created fasta files: ", fasta_count)
    print("Amount of created queries: ", count_added_queries)

    return None


def create_slurmarry(protein_database, input_dir, output_dir):
    """
    Simply creates the slurmarray file for a given input and output directory (each dir with an / at the end!)
    :param protein_database: trivial
    :param input_dir: trivial
    :param output_dir: trivial
    :return: CLCR_slurmarray.slurm file at the cwd
    """

    # creating the .slurm file for the jobarray
    slurm_filename = os.getcwd() + "/" + "CLCR_slurmarray.slurm"
    slurm_file = open(slurm_filename, "w")

    # -d database, -q query, -o output path, -k max hits saved per query seq in output file
    blast_command = "diamond blastx -d " + protein_database + " "\
                    "-q " + input_dir + "temp_in_${SLURM_ARRAY_TASK_ID}.fasta " \
                    "-o " + output_dir + "temp_out_${SLURM_ARRAY_TASK_ID}.txt " \
                    "-k 25 --max-hsps 0 -c 1 -t /dev/shm -F 15 -f 0"

    slurm_file.write("#!/bin/bash\n\n"
                     "#SBATCH --partition=all\n"
                     "#SBATCH --account=praktikant\n"
                     "#SBATCH --cpus-per-task=4\n"
                     "#SBATCH --mem-per-cpu=16gb\n"
                     "#SBATCH --job-name=\"CLCR_run\"\n\n")

    slurm_file.write(blast_command)
    slurm_file.close()

    # sending the job via slurm to the cluster
    #os.system(("sbatch --array=1-" + str(fasta_count) + slurm_filename))

    return None


def count_diamond_hits(project_dir):

    """
    simply count the amount of DIAMOND hits
    :param project_dir: path to the project dir, to find the output dir
    :return: list containing the output regions
    """

    # returns every file in the directory with .out at the end
    output_file_list = glob.glob((project_dir + "/output_dir") + "/temp_out_*.txt")

    output_list = []
    for x in output_file_list:
        current_fasta = open(x, "r")

        for line in current_fasta:
            temp_1 = line.split()
            temp_2 = temp_1[0].split("_")
            output_list.append((int(temp_2[0]), int(temp_2[1]), temp_1[1]))

    #print(output_list)
    print(len(output_list))

    #print("Diamond hits: ", len(output_list))
    # depending on the GUI what to return
    return output_list


def count_all_frameshifts_per_hit(output_dir):

    """
    this version is for the frameshift option, it counts the frameshifts per hit und adds several other information
    :param output_dir: path to the output dir
    :return: list containing the output regions
    """

    # returns every file in the directory with .out at the end
    output_file_list = glob.glob(output_dir + "/temp_out_*.txt")


    # sorting the file path ascending by their file number
    dir_path_len = len(output_dir) + 9  # +9 because temp_out_ consists of 9 chars
    output_file_list = sorted(output_file_list, key=lambda current_path: int(current_path[dir_path_len:-4]))

    output_region_list = []     # containing the Diamond hits with additional information

    # reading out all output files
    for x in output_file_list:
        current_output_file = open(x, "r")

        current_region = ()     # start and end position of the current considered region
        frameshift_count = 0    # counts the amount of frameshifts per query
        hit_detected = False    # becomes true if a hit regions is detected, and becomes False if the current region is
                                # saved in the output_region_list and a new hit is searched
        protein_hit = ""  # containing the protein hit of the current region
        e_value = ""
        bit_score = ""
        similarity_percentage = ""
        query_alignment = ""
        subject_alignment = ""

        # reading out the current output file
        for line in current_output_file:

            # case for a new Query
            if line[0:6] == "Query=":

                # case if the previous region had an Diamond hit
                if hit_detected:

                    query_pos = 0   # for saving the current query position
                    # detect frameshifts in aligned query region
                    for i in query_alignment:
                        # stays True if a putative intron is detected
                        wrong_intron_left = True    # represents the upstream region
                        wrong_intron_right = True   # represents the downstream region
                        if (i == "\\") or (i == "/"):   # frameshift detected

                            # checks if the detected frameshift is a result of a putative Intron
                            # if case for edge positions, to prevent "list out of range errors"
                            if not (((query_pos + 18) > len(query_alignment)) or ((query_pos - 17) < 0)):
                                # used range starts at i-2 to be fault tolerant in intron-exon transition positions
                                for y in range(query_pos-17, query_pos-2):  # upstream region
                                    if subject_alignment[y] != "-":
                                        wrong_intron_left = False

                                for y in range(query_pos+3, query_pos+18):  # downstream region
                                    if subject_alignment[y] != "-":
                                        wrong_intron_right = False

                                # counts the frameshift, if its not caused by a putative intron
                                if not (wrong_intron_right or wrong_intron_left):
                                    frameshift_count += 1

                            # case for frameshifts at the edge of the alignment
                            else:
                                # stays True if the frameshift might be the result of a putative intron
                                wrong_intron = True

                                # case for a downstream edge position
                                if (query_pos + 18) > len(query_alignment):
                                    for y in range(query_pos - 17, query_pos - 2):
                                        if subject_alignment[y] != "-":
                                            wrong_intron = False

                                # case for a upstream edge position
                                else:
                                    for y in range(query_pos + 3, query_pos + 18):
                                        if subject_alignment[y] != "-":
                                            wrong_intron = False

                                if not wrong_intron:
                                    frameshift_count += 1

                        # increasing the value to represent the next position
                        query_pos += 1

                    output_region_list.append([current_region[0], current_region[1], protein_hit, e_value, bit_score,
                                               similarity_percentage, frameshift_count])

                    # reset values
                    frameshift_count = 0
                    query_alignment = ""
                    subject_alignment = ""
                    hit_detected = False

                # setting the new region
                current_line = line.strip().split()[1].split("_")
                current_region = (current_line[0], current_line[1])

            # line containing the protein id
            if line[0] == ">":
                protein_hit = line.split()[0][1:]   # containing the protein hit of the current region
                hit_detected = True

            # line containing score and e-value
            if line[1:6] == "Score":
                current_line = line.strip().split()
                bit_score = current_line[2]
                e_value = current_line[7]

            # line similarity percentage
            if line[1:11] == "Identities":
                current_line = line.strip().split()
                similarity_percentage = current_line[3][1:-2]

            # reading in the aligned sequences
            if line[0:6] == "Query ":
                query_alignment += line.split()[2]
            if line[0:6] == "Sbjct ":
                subject_alignment += line.split()[2]

        # appending the last region if necessary (same as upper part)
        if hit_detected:
            query_pos = 0  # for saving the current query position
            # detect frameshifts in aligned query region
            for i in query_alignment:
                # stays True if a putative intron is detected
                wrong_intron_left = True  # represents the upstream region
                wrong_intron_right = True  # represents the downstream region
                if (i == "\\") or (i == "/"):  # frameshift detected

                    # checks if the detected frameshift is a result of a putative Intron
                    # if case for edge positions, to prevent "list out of range errors"
                    if not (((query_pos + 18) > len(query_alignment)) or ((query_pos - 17) < 0)):
                        # used range starts at i-2 to be fault tolerant in intron-exon transition positions
                        for y in range(query_pos - 17, query_pos - 2):  # upstream region
                            if subject_alignment[y] != "-":
                                wrong_intron_left = False

                        for y in range(query_pos + 2, query_pos + 18):  # downstream region
                            if subject_alignment[y] != "-":
                                wrong_intron_right = False

                        # counts the frameshift, if its not caused by a putative intron
                        if not (wrong_intron_right or wrong_intron_left):
                            frameshift_count += 1

                    # case for frameshifts at the edge of the alignment (always counted)
                    else:
                        # stays True if the frameshift might be the result of a putative intron
                        wrong_intron = True

                        # case for a downstream edge position
                        if (query_pos + 18) > len(query_alignment):
                            for y in range(query_pos - 17, query_pos - 2):
                                if subject_alignment[y] != "-":
                                    wrong_intron = False

                        # case for a upstream edge position
                        else:
                            for y in range(query_pos + 3, query_pos + 18):
                                if subject_alignment[y] != "-":
                                    wrong_intron = False

                        if not wrong_intron:
                            frameshift_count += 1

                # increasing the value to represent the next position
                query_pos += 1

            output_region_list.append([current_region[0], current_region[1], protein_hit, e_value, bit_score,
                                       similarity_percentage, frameshift_count])

        current_output_file.close()

    return output_region_list


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


def read_in_results_3(output_dir, max_detect_dist):
    """
    The function reads in the Diamond blastx results in a given directory. The results for each hit are stored in the
    output_region_list, which is also returned. The contains for every region: Scaffold, start pos., end pos. in
    scaffold, e-value, bitscore, similarity-percentage, list with all detected frameshifts with their positions in the
    query (I = insertion, D = deletion). The second output list (healing_region_list), contains the combined frameshift
    information of all Diamond hits per query, which means the scaffold, the start position in the scaffold, the end
    position in the scaffold, and a list with all positions in the query, where frameshifts are detected and need to
    be healed later . But not all Diamond hits of each Query are included, some are excluded by our overlapping
    heuristic, more detailed information in the documentation and in the code.
    In cases, where the low regions mus be expanded to work as a query. It could happen, that a detected frameshift in
    query is located in the expanded part of the query and not in the original low cov. region. This contradicts the
    concept of the whole program, using low cov. regions that correlate with bad polishing. The max_detect_distance
    defines the distance from a detected frameshift position to the original low cov. region, where a frameshift is
    still considered and not excluded in the further analysis.
    PS: frameshift positions are the python list positions.
    :param output_dir: path to the output directory
    :param max_detect_dist: The max_detect_distance defines the distance from
                            a detected frameshift position to the original low cov. region, where a frameshift is still
                            considered and not excluded in the further analysis. Values around 10 might be reasonable.
    :return: output_region_list, healing_region_list
    """
    # Returns every file in the directory with .out at the end
    output_file_list = glob.glob(output_dir + "/temp_out_*.txt")

    # Sorting the file path ascending by their file number
    dir_path_len = len(output_dir) + 9  # +9 because temp_out_ consists of 9 chars
    output_file_list = sorted(output_file_list, key=lambda current_path: int(current_path[dir_path_len:-4]))

    output_region_list = []  # Containing the Diamond hits with additional information

    # Reading out all output files
    for x in output_file_list:
        current_output_file = open(x, "r")

        # Initialising the parameters
        query_hit_list = []    # Contains all hits of the current query
        current_region = ("initscaffold", 0, 0)  # Scaffold, start and end position of the current considered region
        frameshift_list = []  # Contains the positions of detected frameshifts
        query_alignment = ""
        # The following two parameters are important for the localising of the frameshift in the scaffold, and the later
        # Overlapping heuristic for decide which frameshifts are considered in each query
        query_alignment_start_pos = 0   # Saves the start position of the alignment in the query, initialising it with 0
                                        # is no problem because Diamond starts with counting at 1

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

                    # Appends each protein hit of the query to the output_region_list
                    for current_protein_hit in query_hit_list:

                        # The placed gaps in the query in the alignment are counted for a correct position calculation
                        query_gap_count = 0

                        # Detect frameshifts in aligned query region
                        for query_pos, char in enumerate(current_protein_hit[4]):

                            # Counts how many gaps are inserted into the query before the frameshift position
                            if char == "-":
                                query_gap_count += 1

                            if char == "\\":  # Insertion detected
                                # Function returns True if the frameshift correlates with a putative intron
                                wrong_intron = exclude_putative_transition_frameshift(current_protein_hit[4],
                                                                                      current_protein_hit[5],
                                                                                      query_pos)

                                if not wrong_intron:
                                    # Saves the positions as the python list positions of the first nucleotide in the
                                    # triplet with the frameshift
                                    frameshift_list.append(((current_protein_hit[6] - 1) +
                                                            (3 * (query_pos - query_gap_count)), "I"))

                            elif char == "/":  # Deletion detected
                                # Function returns True if the frameshift correlates with a putative intron
                                wrong_intron = exclude_putative_transition_frameshift(current_protein_hit[4],
                                                                                      current_protein_hit[5],
                                                                                      query_pos)

                                if not wrong_intron:
                                    # Saves the positions as the python list positions of the first nucleotide in the
                                    # triplet with the frameshift
                                    frameshift_list.append(((current_protein_hit[6] - 1) +
                                                            (3 * (query_pos - query_gap_count)), "D"))

                        output_region_list.append([current_region[0], current_region[1], current_region[2],
                                                   current_protein_hit[0], current_protein_hit[1],
                                                   current_protein_hit[2], current_protein_hit[3], frameshift_list,
                                                   current_protein_hit[6], current_protein_hit[7]])

                        # Reset values
                        frameshift_list = []

                    # Reset the the list, the if request is inactive till a new protein hit was found and thus
                    # initialising element is created
                    query_hit_list = []

                # Setting the new region
                current_line = line.strip().split()[1].split("#")
                current_region = (current_line[0], current_line[1], current_line[2])

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

        # Appending the last region if necessary (same as upper part)
        if query_hit_list:
            # Appends the last remaining protein hit of the query
            query_hit_list.append([protein_hit, e_value, bit_score, similarity_percentage, query_alignment,
                                   subject_alignment, query_alignment_start_pos, query_alignment_end_pos])

            # Removes the initialising element of the query
            query_hit_list.pop(0)

            # Appends each protein hit of the query to the output_region_list
            for current_protein_hit in query_hit_list:

                # The placed gaps in the query in the alignment are counted for a correct position calculation
                query_gap_count = 0

                # Detect frameshifts in aligned query region
                for query_pos, char in enumerate(current_protein_hit[4]):

                    # Counts how many gaps are inserted into the query before the frameshift position
                    if char == "-":
                        query_gap_count += 1

                    if char == "\\":  # Insertion detected
                        # Function returns True if the frameshift correlates with a putative intron
                        wrong_intron = exclude_putative_transition_frameshift(current_protein_hit[4],
                                                                              current_protein_hit[5],
                                                                              query_pos)

                        if not wrong_intron:
                            # Saves the positions as the python list positions of the first nucleotide in the
                            # triplet with the frameshift
                            frameshift_list.append(((current_protein_hit[6] - 1) +
                                                    (3 * (query_pos - query_gap_count)), "I"))

                    elif char == "/":  # Deletion detected
                        # Function returns True if the frameshift correlates with a putative intron
                        wrong_intron = exclude_putative_transition_frameshift(current_protein_hit[4],
                                                                              current_protein_hit[5],
                                                                              query_pos)

                        if not wrong_intron:
                            # Saves the positions as the python list positions of the first nucleotide in the
                            # triplet with the frameshift
                            frameshift_list.append(((current_protein_hit[6] - 1) +
                                                    (3 * (query_pos - query_gap_count)), "D"))

                output_region_list.append([current_region[0], current_region[1], current_region[2],
                                           current_protein_hit[0], current_protein_hit[1],
                                           current_protein_hit[2], current_protein_hit[3], frameshift_list,
                                           current_protein_hit[6], current_protein_hit[7]])

                # Reset the list for the next protein hit
                frameshift_list = []

        current_output_file.close()

    """Ich sollte die Art wie die output region list verabeitet/ erstellt wir anpassen. Nicht jeden hit einzeln appenden,
    sondern subliste mit allen hits und den jeweiligen query informationen speichern, also die low cov. region UND die
    query region, welche bei kurzen low cov. abweichen kann. Wenn ich dann die Overlapping heuristic verwende, sollten 
    zwar immer noch die hits nach evalue sortiert einbezogen werden, aber die detectierten frameshifts nur wenn sie nahe
    genug/ in der low cov. region liegen (Query Teil wird von dem berücksichtigten alignment abgedeckt auch wenn die 
    Frameshifts nicht miteinbezogen werden...). Die Anpassung der Frameshift positionen auf die Orginalpositionen und 
    das Aussortieren würde ich schon in dieser funktion machen, somit werden die gefilterten/korrekten positione(welche 
    auch negativ sein können) and die healing funktion übermittelt
    
    Die output region list so neu erstellen, sodass nur hits enthalten sind, welche auch wirklich berücksichtigt wurden,
    sodass die Folgefunktionen die Daten korrekt weiterverarbeiten können."""

    # Determine the frameshift positions in each query, by combining the information of all hits per query.
    # Contains the information of the current query  [scaffold, startpostion_in_scaffold, endposition_in_scaffold,
    # [frameshift position list]]
    current_query = ["#", "#", "#"]
    current_query_region_coverage = []     # Contains the start and end pos. of hit alignments, that
    # healing_region_list contains the final query information for all queries, and functions as output list
    healing_region_list = []

    for diamond_hit in output_region_list:
        # If case for start of new query information, old, now complete, query information could be appended to output
        if diamond_hit[0:3] != current_query[0:3]:
            healing_region_list.append(current_query)       # Appends the final information of the previous region
            # Initialise with new information
            current_query = [diamond_hit[0], diamond_hit[1], diamond_hit[2], copy.copy(diamond_hit[7])]
            current_query_region_coverage = [[int(diamond_hit[8]), int(diamond_hit[9])]]

        # Else case for adding frameshift information of an additional protein hit if possible
        else:
            # Initialise parameters
            region_overlapping = False      # Saves if the new region is overlapping with a already included region
            alignment_start_pos = diamond_hit[8]
            alignment_end_pos = diamond_hit[9]

            # Checks for all already added regions, if the new region is overlapping with one of them
            for prev_region in current_query_region_coverage:
                if ((alignment_end_pos >= prev_region[0]) and (alignment_start_pos <= prev_region[0])) or\
                    ((alignment_end_pos >= prev_region[1]) and (alignment_start_pos <= prev_region[1])) or\
                        ((alignment_start_pos >= prev_region[0]) and (alignment_end_pos <= prev_region[1])):
                    region_overlapping = True       # overlapping region found
                    break

            # Case were no previous region is overlapping with the current region
            if not region_overlapping:
                # Append the new region to the coverage list
                current_query_region_coverage.append([alignment_start_pos, alignment_end_pos])
                """Hier nur appenden, wenn die Frameshifts nahe genug/ in der low cov. region liegen, die 
                Positionsangaben der Frameshifts zusätlich auf den querystart anpassen."""
                current_query[3] += diamond_hit[7]      # add the frameshifts of the new hit to final list

        # Remove the unnecessary start and end positions of the alignment in the output
        diamond_hit.pop(9)
        diamond_hit.pop(8)

    # Append last remaining region
    healing_region_list.append(current_query)
    # Remove the initialising object
    healing_region_list.pop(0)

    return output_region_list, healing_region_list


def heal_assembly_file(healing_region_list, input_fna_path, outut_dir):
    """
    Gets the healing_region_list which contains the frameshift positions in each query, and inserts N's at those
    positions, to heal the reading frame. The healed assembly is afterwards saved in the give output directory.
    The new assembly file doesnt contain comments.
    :param healing_region_list: Contains for each query: scaffold, start pos., end pos. in scaff., frameshift pos. list
    :param outut_dir: directory where the new "healed assembly" is stored (should contain a \ at the last position)
    :param input_fna_path:  File path to the original assembly file
    :return: File path to the new modified assembly file
    """
    # Read in the .fna file for fast access
    input_fna_file = open(input_fna_path)
    scaffold_list = []  # Filled with the sequences of the fna file
    temp_scaffold = []  # Contains the current scaffold
    current_header = ""  # Containing the header of the current scaffold

    # Filling the .fna region list/ reading in the original assembly
    for line in input_fna_file:

        # Sequence line case
        if (line[0] != ">") and (line[0] != ";"):
            for x in line.strip():  # filling up a new sequence
                temp_scaffold.append(x)

        # Sequence header case
        elif line[0] == ">":
            scaffold_list.append([current_header, temp_scaffold])  # appending the previous region
            current_header = line.strip().split()[0][1:]  # reading in the new header
            temp_scaffold = []  # resetting the list/ ready for new region

    # Append last remaining region
    scaffold_list.append([current_header, temp_scaffold])
    input_fna_file.close()
    # Delete the first empty initialising element
    scaffold_list.pop(0)

    print(scaffold_list[0][0], len(scaffold_list[0][1]), scaffold_list[0][1][0:10])

    # Correct the frameshifts
    # Inserting the N's from the end to the beginning, to prevent the necessity of adapting the insertion position in
    # dependency to the previous inserted frameshifts
    # For this, in this first sort, the queries in each scaffold are sorted descending by their start position
    # but the reverse option also sorts the scaffolds descending
    sorted_healing_region_list = sorted(healing_region_list, reverse=True,  key=lambda temp_query: (temp_query[0],
                                                                                                    temp_query[1]))
    print("sorted")
    count = 0

    for query in sorted_healing_region_list:
        # Search the corresponding scaffold
        for scaffold in scaffold_list:
            # Case where matching scaffold is found
            if scaffold[0] == query[0]:
                query_start_pos = int(query[1])
                count += 1
                if count == 1200:
                    print("10% more")
                    count = 0

                # Second sort step, now the frameshift positions in each query are sorted descending
                temp_sorted_query = sorted(query[3], reverse=True,  key=lambda temp_query: temp_query[0])

                for frameshift_position in temp_sorted_query:      # Same usage of reversed() as before
                    # Deletion case
                    if frameshift_position[1] == "D":
                        scaffold[1] = scaffold[1][:(query_start_pos + frameshift_position[0])] + ["N"] + \
                                      scaffold[1][(query_start_pos + frameshift_position[0]):]
                    # Insertion case
                    else:
                        scaffold[1] = scaffold[1][:(query_start_pos + frameshift_position[0])] + ["N"] + ["N"] + \
                                      scaffold[1][(query_start_pos + frameshift_position[0]):]
                break

    # Create the new assembly .fna file path, located in the same dir as the original assembly file
    new_fna_file_path = outut_dir + "new_fna_file.txt"

    print("create new assembly file")

    # Creating the new assembly file (with the inserted N's)
    new_fna_file = open(new_fna_file_path, "w")
    for scaffold in scaffold_list:

        # Write the header line for the scaffold
        new_fna_file.write(">" + scaffold[0] + "\n")

        # Create lines with a length of 80, for a more readable .fasta file
        # Create list with the end positions of each 80 nucleotide long line
        position_list = list(range(80, len(scaffold[1]), 80))
        # Writes the lines in the assembly file
        for list_pos in position_list:
            # The nucleotides are stored as chars in a list and are converted to a string with join
            new_fna_file.write(("".join(scaffold[1][(list_pos - 80):list_pos])) + "\n")

        # Write the last remaining line (range function doesnt include the last elements)
        new_fna_file.write(("".join(scaffold[1][position_list[-1]:len(scaffold[1])])) + "\n")

    new_fna_file.close()

    return new_fna_file_path


def create_gff_adaption_list(healing_region_list):
    """
    Creates a list for the annotation file adaption to the new assembly with the inserted N's, to speed up the process.
    The healing_region_list is used to determine how many frameshifts are inserted till a certain position. Those
    information is saved in the adaption_list, which has a sublist for each scaffold. Each scaffold sublist, contains
    all positions where the amount of inserted N's has changed, which means all positions where frameshifts where
    inserted. For each of this positions is saves how many frameshifts are inserted till this position INCLUDING the
    inserted frameshifts in this position.
    :param healing_region_list: Contains for each query [scaffold, start pos., end pos., [[frameshift position, "I" or
    "D"], [...], ...]
    :return: adaption_list
    """

    # Saves for all scaffolds the positions, when the count of the inserted frameshifts rises. This brings the
    # possibility for a fast check, on how many positions a feature needs to be shifted
    # The List has a sublist for each scaffold, consisting of lists (scaffold_pos, already_inserted_frameshift_count)
    adaption_list = []
    current_scaffold_list = ["####", []]

    # Sort the List after theirs scaffolds and ascending by their start positions in the scaffold
    sorted_healing_region_list = sorted(healing_region_list, key=lambda temp_query: (temp_query[0], temp_query[1]))

    # Parameter initialising
    frameshift_count = 0  # Counts how many frameshifts are inserted in the scaffold, till the current position

    # Fill the list
    for query in sorted_healing_region_list:

        # Case where a new scaffold starts
        if query[0] != current_scaffold_list[0]:
            # Append the old scaffold to the list
            adaption_list.append(current_scaffold_list)

            # Initialise the new scaffold
            frameshift_count = 0
            current_scaffold_list = [query[0], []]

            # Sorts the frameshift positions, of the current query, ascending
            sorted_frameshifts = sorted(query[3])
            # Appends each frameshift to the list
            for single_frameshift in sorted_frameshifts:

                if single_frameshift[1] == "I":     # Insertion case, two additional N's are inserted at this position
                    frameshift_count += 2

                else:       # Deletion case, one additional N is inserted
                    frameshift_count += 1

                # Append the frameshift position and the count of how many frameshifts are inserted till this position
                # to the list, query[1] is the start position of the query in the scaffold and single_frameshift[0] is
                # the position of the frameshift in the in the query
                current_scaffold_list[1].append([(query[1] + single_frameshift[0]), frameshift_count])

        # Case, where the current query is in the same scaffold as the previous region
        else:
            # Sorts the frameshift positions, of the current query, ascending
            sorted_frameshifts = sorted(query[3])
            # Appends each frameshift to the list
            for single_frameshift in sorted_frameshifts:

                if single_frameshift[1] == "I":  # Insertion case, two additional N's are inserted at this position
                    frameshift_count += 2

                else:  # Deletion case, one additional N is inserted
                    frameshift_count += 1

                # Append the frameshift position and the count of how many frameshifts are inserted till this position
                # to the list
                current_scaffold_list[1].append([(query[1] + single_frameshift[0]), frameshift_count])

    # Append the last remaining scaffold
    adaption_list.append(current_scaffold_list)
    # Remove the initialising element
    adaption_list.pop(0)

    return adaption_list


def adapt_gff_file(gff_file_path, adaption_list):
    """
    Adapts a given annotation file to the modified assembly(inserted N's). Shifts the start and end positions of
    annotation features downstream as far a necessary (dependent on how many N's are inserted till those positions).
    The list which is used for the adaption is preprocessed by the create_gff_adaption_list function. Please look at the
    header of this function for a closer description of the used list.
    :param gff_file_path: Path to the .gff file
    :param adaption_list: List which contains sublists for each scaffold with fixed frameshifts. Each scaffold sublist,
    contains all positions where the amount of inserted N's has changed, which means all positions where frameshifts
    where inserted. For each of this positions is saves how many frameshifts are inserted till this position INCLUDING
    the inserted frameshifts in this position.
    :return: None, but creates a new_gff_file.txt with the new annotation
    """

    feature_list = []       # contains all features of the annotation file
    annotation_file = open(gff_file_path)

    for line in annotation_file:
        if line[0] != "#":      # Excludes comments
            line_split = line.split()
            # Corrects the feature start and end positions in the scaffold if necessary
            for scaffold in adaption_list:     # Scaffolds without frameshifts doesn't occur
                if scaffold[0] == line_split[0]:       # Matching scaffold found
                    # Checks how on how many positions the start position needs to be shifted
                    # Searches for the first entry, which is downstream of the start pos., and then uses the frameshift
                    # count of the previous entry

                    # Add initialising objects to fix exception cases
                    temp_frameshift_list = copy.copy(scaffold[1])
                    # Quick solution for the case where no frameshift is downstream of the start or end pos.
                    # Append a last infinite element, which is bigger than all start/end positions
                    temp_frameshift_list.append([float("inf"), 0])
                    # Quick solution for the case where no frameshift is upstream of the start or end pos.
                    # Append a first element, which is smaller than all start/end positions
                    temp_frameshift_list = [[0, 0]] + temp_frameshift_list

                    for list_pos, frameshift in enumerate(temp_frameshift_list):
                        # Case, where the first entry, which is downstream of the start pos., is found
                        if frameshift[0] > int(line_split[3]):

                            # Adds the frameshift count to the start position
                            line_split[3] = str(int(line_split[3]) + temp_frameshift_list[(list_pos - 1)][1])
                            break

                    # Checks how on how many positions the end position needs to be shifted
                    for list_pos, frameshift in enumerate(temp_frameshift_list):
                        # Case, where the first entry, which is downstream of the start pos., is found

                        if frameshift[0] > int(line_split[4]):
                            # Adds the frameshift count to the start position
                            line_split[4] = str(int(line_split[4]) + temp_frameshift_list[(list_pos - 1)][1])
                            break

                    break

            feature_list.append("\t".join(line_split))

    annotation_file.close()

    # Create the new annotation file
    new_gff_file_path = gff_file_path[:(gff_file_path.rfind("/") + 1)] + "new_gff_file.gff"

    new_gff_file = open(new_gff_file_path, "w")

    # Write the adapted features in the new file
    for feature_line in feature_list:
        new_gff_file.write(feature_line + "\n")

    new_gff_file.close()

    return None


def add_gff_information(gff_file, output_region_list, fna_file):

    """
    adding to each detected region the overlapping feature information of the annotation file
    :param gff_file: path to gff file
    :param output_region_list: list with the output regions
    :param fna_file: path to the fna file
    :return: the output_region_list with the additional feature information
    """

    # reading out the different .fasta/.fna sequences, to assign th correct positions in the complete sequence to the
    # .gff regions, for a correct assignment from low cov. region to .gff region

    # initialising the parameters
    fna_region_list = []
    current_fna_region = "#"
    current_position = 0
    region_start = 0
    input_file_2 = open(fna_file)

    print("reading in .fna file")
    # reading in the start and the end of each region
    for line in input_file_2:

        # detecting new region, appending the old region to the list
        """ACHTUNG AUSLESEN DER SCAFFOLDS SEHR SPEZIFISCH (HEADER ZEILE KANN ANDERES FORMAT HABEN!!!!!)"""
        if line[0] == ">":
            fna_region_list.append((current_fna_region, region_start, current_position))
            region_start = current_position
            current_fna_region = line[1:(len(line))].strip()

        else:
            if line[0] != ";":
                current_position += len(line.strip())

    # appending the last remaining element
    fna_region_list.append((current_fna_region, region_start, current_position))
    # deleting the initialising element
    fna_region_list.pop(0)

    # read the regions out of the gff file and saving them in region_list
    input_file = open(gff_file)
    gff_region_list = []  # containing tuple like (region_start, region_end, region_type)
    # counting the regions where no parent seq. was found at the .fna file (there should be a parent regions, and if not
    # something went wrong!!!)
    parent_not_found = 0

    print("reading in .gff file")
    # reading the exon regions out of the gff
    for line in input_file:
        if line[0] != "#":
            # for excluding the false line (lines without gff table)
            if len(line.split()) > 7:
                if (line.split()[2] == "exon") or (line.split()[2] == "intron"):

                    # calculating the positions in the whole sequence, using the scaffold region positions
                    parent_sequence = line.split()[0]

                    # searching the corresponding parent sequence for each .gff region
                    parent_found = False    # value containing the information if a parent sequence was found (there
                    # should be a parent seq. if not there ist some thing wrong with the .gff)
                    for parent in fna_region_list:
                        if parent[0] == parent_sequence:
                            parent_reg_start = parent[1]
                            gff_reg_start = parent_reg_start + int(line.split()[3])
                            gff_reg_end = parent_reg_start + int(line.split()[4])
                            parent_found = True
                            break

                    # appending the .gff region start and end in the complete sequence
                    if parent_found:
                        gff_region_list.append((gff_reg_start, gff_reg_end, line.split()[2]))

                    # counting the problem regions where no parent seq. was found
                    else:
                        parent_not_found += 1

    input_file.close()

    print("lenght of gff_region_list: ", len(gff_region_list))

    # sorting the unsorted region_list ascending by the regions end positions
    gff_region_list = sorted(gff_region_list, key=lambda current_tuple: current_tuple[1])
    # invert the sorted list to sort them descending by their region end positions
    gff_region_list.reverse()       # no runtime differences to instead using the ascending list in a for loop

    print("combining the information")
    # assign regions of output_region_list to regions in region_list if possible
    for current_out_region in output_region_list:

        current_gff_region = 0
        out_reg_start = int(current_out_region[0])      # starting position of the current "output" region
        out_reg_end = int(current_out_region[1])       # end position of the current "output" region
        region_type = []  # for saving the output region type, if there is a match with a gff region

        # runs as long as the start of the output region is smaller than the end of the current gff_region
        while out_reg_start <= gff_region_list[current_gff_region][1]:
            gff_reg_start = gff_region_list[current_gff_region][0]
            gff_reg_end = gff_region_list[current_gff_region][1]

            # checks if either the output regions start ot the end lays in a gff region, or the gff region is
            # completely covered
            if ((out_reg_end >= gff_reg_start) and (out_reg_end <= gff_reg_end)) or \
                ((out_reg_start >= gff_reg_start) and (out_reg_start <= gff_reg_end)) or \
                    ((out_reg_start <= gff_reg_start) and (out_reg_end >= gff_reg_end)):
                # case if region type is empty
                gff_type = gff_region_list[current_gff_region][2]   # current gff region type
                if not region_type:
                    # appends the gff region type
                    region_type.append(gff_type)

                # case if the output region lays in more than one gff region
                else:
                    # checks if the region is already assigned to that gff region type
                    new_type = True
                    for x in region_type:
                        if x == gff_type:
                            new_type = False

                    if new_type:    # case if a new type was found
                        region_type.append(gff_type)

            current_gff_region += 1
            # break condition
            if current_gff_region == len(gff_region_list):
                break

        # case if the output region was assigned to a gff_region
        if region_type:
            # case if both region types are covered
            if len(region_type) == 2:
                current_out_region.append("exon_intron")
            # case if only one region type is covered
            else:
                current_out_region.append(region_type[0])

        # case if there was no assignment, so the region is noncoding
        else:
            current_out_region.append("noncoding")

    print("Regions with no parentsequences (This shouldn't happen!): ", parent_not_found)

    return output_region_list


def add_gff_information_mod(gff_file, output_region_list, fna_file):

    """
    function when no intron regions are included in the -gff file (instead using gene regions)
    adding to each detected region the overlapping feature information of the annotation file
    :param gff_file: path to gff file
    :param output_region_list: list with the output regions
    :param fna_file: path to the fna file
    :return: the output_region_list with the additional feature information
    """

    # reading out the different .fasta/.fna sequences, to assign th correct positions in the complete sequence to the
    # .gff regions, for a correct assignment from low cov. region to .gff region

    # initialising the parameters
    fna_region_list = []
    current_fna_region = "#"
    current_position = 0
    region_start = 0
    input_file_2 = open(fna_file)

    print("reading in .fna file")
    # reading in the start and the end of each region
    for line in input_file_2:

        # detecting new region, appending the old region to the list
        """ACHTUNG AUSLESEN DER SCAFFOLDS SEHR SPEZIFISCH (HEADER ZEILE KANN ANDERES FORMAT HABEN!!!!!)"""
        if line[0] == ">":
            fna_region_list.append((current_fna_region, region_start, current_position))
            region_start = current_position
            current_fna_region = line[1:(len(line))].strip()

        else:
            if line[0] != ";":
                current_position += len(line.strip())

    input_file_2.close()

    # appending the last remaining element
    fna_region_list.append((current_fna_region, region_start, current_position))

    # deleting the initialising element
    fna_region_list.pop(0)

    # read the regions out of the gff file and saving them in gff_region_list
    input_file = open(gff_file)
    gff_region_list = []  # containing tuple like (region_start, region_end, region_type)
    # counting the regions where no parent seq. was found at the .fna file (there should be a parent regions, and if not
    # something went wrong!!!)
    parent_not_found = 0

    print("reading in .gff file")
    # reading the exon regions out of the gff
    for line in input_file:
        if line[0] != "#":
            # for excluding the false line (lines without gff table)
            if len(line.split()) > 7:
                if (line.split()[2] == "exon") or (line.split()[2] == "gene"):

                    # calculating the positions in the whole sequence, using the scaffold region positions
                    parent_sequence = line.split()[0]

                    # searching the corresponding parent sequence for each .gff region
                    parent_found = False    # value containing the information if a parent sequence was found (there
                    # should be a parent seq. if not there ist some thing wrong with the .gff)
                    for parent in fna_region_list:
                        if parent[0] == parent_sequence:
                            parent_reg_start = parent[1]
                            gff_reg_start = parent_reg_start + int(line.split()[3])
                            gff_reg_end = parent_reg_start + int(line.split()[4])
                            parent_found = True
                            break

                    # appending the .gff region start and end in the complete sequence
                    if parent_found:

                        # CDS is changed to intron, to fit for the later assignment approach (later no type
                        # change needed)
                        if line.split()[2] == "gene":
                            gff_region_list.append((gff_reg_start, gff_reg_end, "intron/UTR"))

                        # normal exon case
                        else:
                            gff_region_list.append((gff_reg_start, gff_reg_end, line.split()[2]))

                    # counting the problem regions where no parent seq. was found
                    else:
                        parent_not_found += 1

    input_file.close()

    print("lenght of gff_region_list: ", len(gff_region_list))

    for x in gff_region_list:
        if (x[0] == 5245342) or (x[0] == 5245414) or (x[0] == 5245710):
            print(x)

    # sorting the unsorted region_list ascending by the regions end positions
    gff_region_list = sorted(gff_region_list, key=lambda current_tuple: current_tuple[1])
    # invert the sorted list to sort them descending by their region end positions
    gff_region_list.reverse()       # no runtime differences to instead using the ascending list in a for loop

    print("assingning annotation features if possible")
    # assign regions of output_region_list to regions in region_list if possible

    progress_count = 0
    len_out_reg_4 = int(len(output_region_list) / 4)

    for current_out_region in output_region_list:

        # printing the progress
        if progress_count == len_out_reg_4:
            print("25% reached")
        if progress_count == (len_out_reg_4 * 2):
            print("50% reached")
        if progress_count == (len_out_reg_4 * 3):
            print("75% reached")


        current_gff_region = 0
        out_reg_start = int(current_out_region[0])      # starting position of the current "output" region
        out_reg_end = int(current_out_region[1])       # end position of the current "output" region
        region_type = []  # for saving the output region type, if there is a match with a gff region

        # runs as long as the start of the output region is smaller than the end of the current gff_region
        while out_reg_start <= gff_region_list[current_gff_region][1]:
            gff_reg_start = gff_region_list[current_gff_region][0]
            gff_reg_end = gff_region_list[current_gff_region][1]

            # checks if either the output regions start or the end lays in a gff region, or the gff region is
            # completely covered
            if ((out_reg_end >= gff_reg_start) and (out_reg_end <= gff_reg_end)) or \
                ((out_reg_start >= gff_reg_start) and (out_reg_start <= gff_reg_end)) or \
                    ((out_reg_start <= gff_reg_start) and (out_reg_end >= gff_reg_end)):
                # case if region type is empty
                gff_type = gff_region_list[current_gff_region][2]   # current gff region type
                if not region_type:
                    # appends the gff region type
                    region_type.append((gff_type, gff_reg_start, gff_reg_end))

                # case if the output region lays in more than one gff region
                else:
                    # checks if the region is already assigned to that gff region type
                    new_type = True
                    for x in region_type:
                        if x[0] == gff_type:
                            new_type = False

                    if new_type:    # case if a new type was found
                        region_type.append((gff_type, gff_reg_start, gff_reg_end))

            current_gff_region += 1
            # break condition
            if current_gff_region == len(gff_region_list):
                break

        # case if the output region was assigned to a gff_region
        if region_type:
            # case if the region covers exon and cds
            if len(region_type) == 2:

                # assigning the gff type region lenght
                for reg_type in region_type:
                    if reg_type[0] == "exon":
                        exon_start = reg_type[1]
                        exon_end = reg_type[2]
                        break

                # case where the low cow. regions lay completely in the exon
                if (out_reg_start >= exon_start) and (out_reg_end <= exon_end):
                    current_out_region.append("exon")

                # case where an exon and an intron is covered by the low cov. region
                else:
                    current_out_region.append("exon_intron")

            # case if only one region type is covered
            else:
                current_out_region.append(region_type[0][0])

        # case if there was no assignment, so the region is noncoding
        else:
            current_out_region.append("noncoding")

        progress_count += 1

    print("Regions with no parentsequences (This shouldn't happen!): ", parent_not_found)
    print("lenght output_region_list: ", len(output_region_list))

    # adding the distance to the next gene
    cds_region_list = []  # list containing all CDS Regions out of the .gff file

    # appending all CDS regions in the new list
    for region in gff_region_list:
        if region[2] == "intron":   # remember, CDS was changed to intron
            cds_region_list.append([region[0], region[1], "CDS"])

    print("adding distance to next CDS information...")
    progress_count = 0
    len_out_reg_4 = int(len(output_region_list)/4)


    # adding the information to every output region
    for current_region in output_region_list:

        if progress_count == len_out_reg_4:
            print("25% reached")
        if progress_count == (len_out_reg_4 * 2):
            print("50% reached")
        if progress_count == (len_out_reg_4 * 3):
            print("75% reached")


        shortest_distance = cds_region_list[0][1]  # initialising the variable with the largest possible distance
        current_cds_region = 0      # count for saving the currently used cds region (tuple position in list)
        out_reg_start = int(current_region[0])  # starting position of the current "output" region

        # the regions in the cds_region_list are still sorted descending by their end position
        # runs as long as the start of the output region is smaller than the end of the current cds_region
        while out_reg_start < cds_region_list[current_cds_region][1]:
            # case that the current_region is upstream from the cds_region, AND not overlapping
            if int(current_region[1]) < cds_region_list[current_cds_region][0]:
                temp_distance = cds_region_list[current_cds_region][0] - int(current_region[1])  # distance between the considere regions
                # case that the new distance from the current_region end to the cds_region start is smaller then the
                # current one
                if temp_distance < shortest_distance:
                    shortest_distance = temp_distance
            current_cds_region += 1

            if current_cds_region == len(cds_region_list):
                break

        # consideres the cds regions with a start upstream of the output region start, the first of those regions,
        # which is not overlapping with the output region, is potentially closer than the closest downstream region
        cds_region_list_len = len(cds_region_list)
        while True:
            if current_cds_region == cds_region_list_len:
                break
            # case that the first upstream cds region, which is not overlapping with the output region is found
            if cds_region_list[current_cds_region][1] < out_reg_start:
                temp_distance = int(current_region[0]) - cds_region_list[current_cds_region][1]
                if temp_distance < shortest_distance:
                    shortest_distance = temp_distance
                break
            current_cds_region += 1

        # append the calculated distance
        current_region.append(shortest_distance)

        progress_count += 1

    return output_region_list


def create_combined_output_file(output_region_list, fna_file):

    """
    Creating the output .txt file, either from regions without or with gff information
    reading out the different .fasta/.fna sequences, to transform the complete sequence positions to the position in
    the specific scaffolds
    :param output_region_list: list containing the output regions
    :param fna_file: path to the fna_file
    :return: nothing
    """

    # initialising the parameters
    fna_region_list = []
    current_fna_region = "#"
    current_position = 0
    region_start = 0
    input_file_2 = open(fna_file)

    print("reading in .fna file...")
    # reading in the start and the end of each region
    for line in input_file_2:

        # detecting new region, appending the old region to the list
        """ACHTUNG AUSLESEN DER SCAFFOLDS SEHR SPEZIFISCH (HEADER ZEILE KANN ANDERES FORMAT HABEN!!!!!)"""
        if line[0] == ">":
            fna_region_list.append((current_fna_region, region_start, current_position))
            region_start = current_position
            current_fna_region = line[1:(len(line))].strip()

        else:
            if line[0] != ";":
                current_position += len(line.strip())

    # appending the last remaining element
    fna_region_list.append((current_fna_region, region_start, current_position))
    # deleting the initialising element
    fna_region_list.pop(0)

    # transforming the complete sequence positions to the position in the specific scaffolds
    new_output_region_list = []
    for output_region in output_region_list:
        output_region_start = int(output_region[0])
        output_region_end = int(output_region[1])
        output_region_scaffold = ""     # assigned scaffold of the current output region
        # searching the fitting  scaffold:
        for scaffolds in fna_region_list:
            # case if a fitting scaffold was found
            if (output_region_start >= scaffolds[1]) and (output_region_start <= scaffolds[2]):
                # modifying the values
                # assigned scaffold of the current output region
                output_region_scaffold = scaffolds[0]
                # start position in the assigned scaffold of the current output region
                output_region_start = output_region_start - scaffolds[1]
                # end position in the assigned scaffold of the current output region
                output_region_end = output_region_end - scaffolds[1]

        new_output_region_list.append([output_region_scaffold, output_region_start, output_region_end]
                                      + output_region[2:])

    # creating the output file
    output_file_path = os.getcwd() + "/output_file.txt"

    output_file = open(output_file_path, "w")

    for region in new_output_region_list:
        new_line = ""
        for region_info in region:
            new_line += (str(region_info) + "\t")
        output_file.write((new_line + "\n"))

    output_file.close()


def filter_output_file(output_file_path, identity_threshold, frameshift_filter):

    """
    filtering the output by throwing out the regions below a given identity threshold and if configured excluding regions
    without a frameshift
    :param output_file_path: output_file_path
    :param identity_threshold: identity threshold
    :param frameshift_filter: if TRUE regions without frameshifts are excluded
    :return: nothing
    """

    # containing the regions with an value higher than the threshold
    filtered_output_regions = []

    output_file = open(output_file_path)

    if frameshift_filter:
        for line in output_file:
            if (int(line.split()[6][:-1]) >= identity_threshold) and (int(line.split()[7]) > 0):
                filtered_output_regions.append(line)
    else:
        for line in output_file:
            if int(line.split()[6][:-1]) >= identity_threshold:
                filtered_output_regions.append(line)

    output_file.close()

    print("regions after filtering with threshold ", identity_threshold, ": ", len(filtered_output_regions))

    # creating the output file
    new_file_path = os.getcwd() + "/filtered_output_file.txt"
    new_file = open(new_file_path, "w")

    for line in filtered_output_regions:
        new_file.write(line)

    new_file.close()

    return None


def transform_txt_to_gff(txt_file_path):

    """
    creating a gff file out of the original txt output file
    :param txt_file_path: file path
    :return: nothing
    """

    txt_file = open(txt_file_path)

    # reading in the txt regions
    gff_line_list = []

    """for line in txt_file:

        line_t = line.strip().split()       # temporary list
        modified_line = line_t[0] + "\t" + "." + "\t" + "LowCov" + "\t" + line_t[1] + "\t" + line_t[2] + "\t" \
                        + ".\t.\t.\t" + "protein_hit \"" + line_t[3] + "\"; evalue \"" + line_t[4] + "\"; bitscore \"" \
                        + line_t[5] + "\"; identities \"" + line_t[6] + "\"; frameshift_count \"" + line_t[7] \
                        + "\"; assigned_gff_region_type \"" + line_t[8] + "\"; distance_to_next_CDS \"" + line_t[9] \
                        + "\";\n"

        gff_line_list.append(modified_line)"""

    for line in txt_file:

        line_t = line.strip().split()       # temporary list
        modified_line = line_t[0] + "\t" + "." + "\t" + "LowCov" + "\t" + line_t[1] + "\t" + line_t[2] + "\t" \
                        + ".\t.\t.\t" + "protein_hit=\"" + line_t[3] + "\";evalue=\"" + line_t[4] + "\";bitscore=\"" \
                        + line_t[5] + "\";identities=\"" + line_t[6] + "\";frameshift_count=\"" + line_t[7] \
                        + "\";assigned_gff_region_type=\"" + line_t[8] + "\";distance_to_next_CDS=\"" + line_t[9] \
                        + "\";\n"

        gff_line_list.append(modified_line)

    txt_file.close()

    # creating the gff file
    new_gff_file_path = os.getcwd() + "/output_gff_file.gff"
    new_gff_file = open(new_gff_file_path, "w")

    for line in gff_line_list:
        new_gff_file.write(line)

    new_gff_file.close()

    return None


def main():
    print("database comparison main executed")


if __name__ == '__main__':
    main()
