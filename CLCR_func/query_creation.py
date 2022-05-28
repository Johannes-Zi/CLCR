#!python3

"""This file contains the functions for the detection of low coverage regions in the pbc file of an assembly"""
__author__ = "6947325: Johannes Zieres"
__credits__ = ""
__email__ = "johannes.zieres@gmail.com"

import time
import os
import argparse


def detect_regions(cov_file_path, cov_start, cov_end, verbose_func):
    """
    Detecting the low coverage regions in the per base coverage file. Saves for each low coverage region the scaffold
    and the start and end position in the scaffold. Low cov. regions that are spanning two or more scaffolds are
    excluded. Returns the low cov regions as python list positions(-1). A low cov. region is detected when the cov. is
    equal or lower than cov_start. And a low cov. region is ended when the coverage is higher than cov_end.
    :param cov_file_path: path of input coverage file
    :param cov_start: threshold for detecting a low coverage area
    :param cov_end: threshold for leaving a low coverage area
    :param verbose_func: verbose option
    :return: list with the low coverage regions
    """

    import_file = open(cov_file_path)  # Opens the coverage file
    current_position = 0  # Saves the current position
    low_cov = False  # Saves the current coverage situation
    low_cov_regions = []  # Low coverage regions List with the current_scaffold list of each scaffold
    region_start = int  # Saves the start point for a low coverage region

    # Saves the scaffold of the first low cov position of the current low coverage region
    start_scaffold = str    # Parameter is for excluding low cov. regions that are spanning multiple scaffolds
    current_scaffold = ["scaffoldname", [(0, 10)]]    # Saves the name of the scaffold and all low cov. regions in them.

    # Detects the low coverage regions and saves them in low_coverage_regions list
    for line in import_file:
        splitted_line = line.split()

        # New scaffold starts, old scaffold information is appended to the final output list
        if splitted_line[0] != current_scaffold[0]:

            if verbose_func:
                print("\033[A                             \033[A", flush=True)
                print("Current scaffold: ", splitted_line[0], flush=True)

            # Appends the old scaffold to the final output list if low cov. reg. are found in the finished scaffold
            if current_scaffold[1]:
                low_cov_regions.append(current_scaffold)
            current_scaffold = [splitted_line[0], []]   # Resets the list and initialise the List for the new scaffold
            current_position = 0    # Resets the current position to 0, because a new scaffold starts

        if low_cov:  # Case for being in a low coverage region
            if int(splitted_line[2]) > cov_end:     # Checks if the low cov. region ends
                low_cov = False

                # Checks that the start and end position of the low cov. regions lays in the same scaffold
                if current_scaffold[0] == start_scaffold:
                    current_scaffold[1].append((region_start, (current_position - 1)))

        else:  # Case for being in a good coverage region
            if int(splitted_line[2]) <= cov_start:      # Checks if a new low cov. region starts
                low_cov = True
                region_start = current_position     # Saves the start position of the current low cov. region
                start_scaffold = splitted_line[0]   # Saves the scaffold of the start position

        current_position += 1

    import_file.close()  # Closes the input file

    # Appends the last remaining scaffold if low coverage regions where detected in them
    if current_scaffold[1]:
        low_cov_regions.append(current_scaffold)

    low_cov_regions.pop(0)      # Removes the initialising object

    return low_cov_regions


def create_low_cov_tsv_file(low_cov_regions, output_file_path):
    """
    This function creates a output tsv file with all originally detected low coverage functions before the merging of
    close  together low cov. regions, and the step where short low cov. regions are expanded.
    This is done for a later comparison of the by Diamond detected frameshifts to the original low cov. regions, thus
    only frameshift positions in low cov. regions are considered in the healing process.
    The low cov regions are stored in the format: scaffold \t start-position-in-scaffold \t end-pos-in-scaff + \n.
    :param low_cov_regions: output of the detect_regions() function
    :param output_file_path: path where the low cov. regions storage file should be saved (including the file name!)
    :return: only the output file
    """

    # Create new low cov storage file
    low_cov_storage_file = open(output_file_path, "w")

    # Store the low cov data in the file
    for scaffold in low_cov_regions:

        current_scaffold_name = scaffold[0]

        # Save all low cov. reg. in separate lines
        for low_cov_region in scaffold[1]:
            low_cov_storage_file.write(current_scaffold_name + "\t" + str(low_cov_region[0]) + "\t" +
                                       str(low_cov_region[1]) + "\n")

    low_cov_storage_file.close()

    return None


def read_in_low_cov_tsv_file(storage_tsv_file_path):
    """
    This functions simply reads in the .tsv file, which was created by the create_low_cov_tsv_file() function.
    The tsv stores all original low cov. regions of an current run.
    The function returns the low_cov_regions list that was used to create the .tsv file. (same as the output of the
    detect_regions function).
    For more information see the .tsv creation function.
    :param storage_tsv_file_path: path to the .tsv file (including the file name!)
    :return: low_cov_regions
    """

    low_cov_regions = []
    current_scaff_name = "initialising_scaffold"
    current_scaff_low_cov_reg_list = []

    # Open low cov storage file
    low_cov_storage_file = open(storage_tsv_file_path, "r")

    for line in low_cov_storage_file:
        line_split = line.strip().split()
        scaff_name = line_split[0]
        low_cov_start_pos = int(line_split[1])
        low_cov_end_pos = int(line_split[2])

        # Case for a new scaffold
        if scaff_name != current_scaff_name:

            # Append previous scaffold to the list
            low_cov_regions.append([current_scaff_name, current_scaff_low_cov_reg_list])

            # Initialise new scaffold lists
            current_scaff_name = scaff_name
            current_scaff_low_cov_reg_list = [(low_cov_start_pos, low_cov_end_pos)]

        # Case for low cov region in the same scaffold as before
        else:
            # Append low cov region to the list
            current_scaff_low_cov_reg_list.append((low_cov_start_pos, low_cov_end_pos))

    # Append the last remaining scaffold
    low_cov_regions.append([current_scaff_name, current_scaff_low_cov_reg_list])

    # Remove initialising scaffold of the list
    low_cov_regions.pop(0)

    return low_cov_regions


def merge_close_regions(input_scaffold_list, merge_distance):

    """
    Receives the output list of the detect regions function and merges all low cov. regions in each scaffold together if
    the are closer together than the given merge_distance parameter.
    :param input_scaffold_list: input list, which is containing the scaffold lists with the low coverage regions
    :param merge_distance: distance threshold, for merging regions which are closer together than the threshold
    :return: region list with the merged regions
    """

    output_list = []    # Contains the modified scaffold lists with low cov. regions of each scaffold

    # Iterates trough the scaffolds
    for scaffold_list in input_scaffold_list:

        new_scaffold_list = [scaffold_list[0], []]      # Saves the new merged low cov. regions

        # Initialising the parameters with the first low cov region
        current_start = scaffold_list[1][0][0]
        current_end = scaffold_list[1][0][1]

        # Merging the regions of the current scaffold
        for low_cov_reg in scaffold_list[1]:

            # Case in which the start position of the current region ist closer to the end position of the previous
            # region than the merge distance, so those two regions must be merged
            if (low_cov_reg[0] - merge_distance) <= current_end:
                current_end = low_cov_reg[1]

            # Case where the current region don't have to be merged with the previous region, which means no following
            # region has to be merged with the previous either. So the previous region is complete and could be appended
            # to the scaffold list
            else:
                new_scaffold_list[1].append((current_start, current_end))

                # Reset the parameters and initialise them with the new low cov. region
                current_start = low_cov_reg[0]
                current_end = low_cov_reg[1]

        # Append the last region
        new_scaffold_list[1].append((current_start, current_end))

        # Append the final scaffold list to the output list
        output_list.append(new_scaffold_list)

    return output_list


def query_files_creation(list_with_scaffold_specific_low_cov_reg_lists, fna_file_path, min_length, query_dir,
                         seq_per_fasta):

    """
    This function creates the .fasta query files for a diamond slurm job, out of a input list with the low coverage
    regions of each and a path to the matching fna file to read out the sequences of the low coverage regions.
    The saved positions of the queries/low cov. regions in each scaffold are the normal positions, not python list
    positions(-1) !!!
    The query header consists of:
    >the scaffold # the start position of the merged low cov region in scaffold # end pos of the merged.. # start pos of
     the query in scaff. # end pos of the query..
    :param list_with_scaffold_specific_low_cov_reg_lists: list consisting of scaffold specific sublists with tuples
                                                          containing the start and he end of a region for database
                                                          comparison (output of detect_regions, merge_regions or
                                                          combine_regions_multiple_scaffolds).
    :param fna_file_path: path of the fna file for receiving the base sequence
    :param min_length: expanding the detected regions to 500 if they are smaller, for a better working on diamond
    :param query_dir: directory where the query files are stored in
    :param seq_per_fasta: max amount of queries stored per fasta, The fisrt and the last file could contain in special
                          cases more queries
    :return fasta_count (amount of created query files), count_added_queries (amount of created queries)
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
    fasta_count = 1     # for individual naming at creating of the .fasta files
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

                    # Creating .fasta file, if the current file is full (defined by seq_per_fasta)
                    if already_added_queries >= seq_per_fasta:
                        fasta_count += 1
                        current_fasta.close()
                        already_added_queries = 0
                        new_fasta = query_dir + "/" + "temp_in_" + str(fasta_count) + ".fasta"
                        current_fasta = open(new_fasta, "w")

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

                    # Writing the query header in the fasta file
                    current_fasta.write((">" + scaffold_name + "#" + str(current_low_cov_list[current_tuple][0]) +
                                         "#" + str(current_low_cov_list[current_tuple][1]) +
                                         "#" + str(region_start) + "#" + str(region_end) + "\n"))

                    # Writing the query sequence in the fasta file
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

                    # Creating .fasta file, if the current file is full (defined by seq_per_fasta)
                    if already_added_queries >= seq_per_fasta:
                        fasta_count += 1
                        current_fasta.close()
                        already_added_queries = 0
                        new_fasta = query_dir + "/" + "temp_in_" + str(fasta_count) + ".fasta"
                        current_fasta = open(new_fasta, "w")

                    # Calculation the length a region must be expanded
                    region_length = current_low_cov_list[current_tuple][1] - current_low_cov_list[current_tuple][0]

                    if region_length < min_length:
                        expand_len = int((min_length - region_length) / 2)
                    else:
                        expand_len = 0

                    # Writes the query header
                    current_fasta.write((">" + scaffold_name + "#" + str(current_low_cov_list[current_tuple][0]) +
                                         "#" + str(current_low_cov_list[current_tuple][1]) +
                                         "#" + str(current_low_cov_list[current_tuple][0] - expand_len) +
                                         "#" + str(current_low_cov_list[current_tuple][1] + expand_len) + "\n"))

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

                # Possible last cases where the region cant be expanded 250 to the right,
                # because (the endposition(of the region) + 250) > sequence length
                while current_tuple < len(current_low_cov_list):

                    # Creating .fasta file, if the current file is full (defined by seq_per_fasta)
                    if already_added_queries >= seq_per_fasta:
                        fasta_count += 1
                        current_fasta.close()
                        already_added_queries = 0
                        new_fasta = query_dir + "/" + "temp_in_" + str(fasta_count) + ".fasta"
                        current_fasta = open(new_fasta, "w")

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

                    #  Writing the query header in the fasta file
                    current_fasta.write((">" + scaffold_name + "#" + str(current_low_cov_list[current_tuple][0]) +
                                         "#" + str(current_low_cov_list[current_tuple][1]) +
                                         "#" + str(region_start) + "#" + str(region_end) + "\n"))
                    # Writing the query sequence in the fasta file
                    current_fasta.write(("".join(current_scaffold[1][region_start:region_end]) + "\n"))
                    current_tuple += 1
                    count_added_queries += 1
                    already_added_queries += 1

                break       # Break because there is only one matching scaffold sequence for each scaffold

    # Close the last .fasta file
    current_fasta.close()

    # User information
    # print("Amount of created fasta files: ", fasta_count)
    # print("Amount of created queries: ", count_added_queries)

    return fasta_count, count_added_queries


def create_queries(args):
    """
    Function for the query creation part, including low coverage region detection, merging of low cov. regions that are
    close together and the query file creation.
    Creates a new storage_dir if necessary, creates a input and output dir (overwrite if already present)
    :param args: Arguments from parser

    The following parameters are handed over in the form of args
    -> project_dir: Should include / at the end
    -> cov_file_path: path to the coverage file
    -> assembly_file: path to the assembly file
    -> low_cov_start: threshold for introducing a low cov region
    -> low_cov_end: threshold for ending a low cov regions
    -> min_query_len: minimum query length in bp
    -> queries_per_file: amount of queries placed per diamond input file

    :return: None
    """

    # Initialise args as values for a better overview
    project_dir = args.project_dir
    cov_file_path = args.cov_file_path
    assembly_file = args.assembly_file
    low_cov_start = args.low_cov_start
    low_cov_end = args.low_cov_end
    min_query_len = args.min_query_len
    queries_per_file = args.queries_per_file
    verbose_func = args.verbose

    if verbose_func:
        print("#### Query creation function called! ####", flush=True)

    # Stores the relevant data of the current run
    run_information = ["CLCR create_queries run \t\t" + time.ctime(time.time())]      # Initialising

    start_time = time.time()

    if verbose_func:
        print("### Detect low cov regions ###", flush=True)
    # Detect the low coverage regions
    low_cov_regions = detect_regions(cov_file_path, low_cov_start, low_cov_end, verbose_func)

    # Count the amount of low cov. regions before merging
    region_count = 0
    for scaffold in low_cov_regions:
        region_count += len(scaffold[1])
    run_information.append("Detected low cov. regions:\t\t" + str(region_count))

    storage_files_dir_path = project_dir + "storage_files/"
    # Create a storage dir, if it not exists
    os_command = " if ! [ -d " + storage_files_dir_path + " ] ; then mkdir " + storage_files_dir_path + "; fi"
    os.system(os_command)

    if verbose_func:
        print("### Create original_low_cov_regions.tsv ###", flush=True)
    # Create the "original" low cov regions file
    low_cov_storage_tsv_path = storage_files_dir_path + "original_low_cov_regions.tsv"
    create_low_cov_tsv_file(low_cov_regions, low_cov_storage_tsv_path)

    if verbose_func:
        print("### Merge low cov regions ###", flush=True)
    # Merge low cov regions, that are closer together than the min_query_length/2
    low_cov_regions = merge_close_regions(low_cov_regions, int(min_query_len/2))

    # Count the regions after merging
    region_count = 0
    for scaffold in low_cov_regions:
        region_count += len(scaffold[1])
    run_information.append("Low cov. regions after merging:\t" + str(region_count))

    # Delete the query_files directory with all files in it, if it already exists
    query_files_dir_path = project_dir + "query_files/"
    os_command = " if [ -d " + query_files_dir_path + " ] ; then rm -r " + query_files_dir_path + "; fi"
    os.system(os_command)

    # Create the new query_files directory
    os_command = "mkdir " + query_files_dir_path
    os.system(os_command)

    # Delete the diamond_output directory with all files in it, if it already exists
    diamond_output_dir_path = project_dir + "diamond_output/"
    os_command = " if [ -d " + diamond_output_dir_path + " ] ; then rm -r " + diamond_output_dir_path + "; fi"
    os.system(os_command)

    # Create the new diamond_output directory
    os_command = "mkdir " + diamond_output_dir_path
    os.system(os_command)

    if verbose_func:
        print("### Create new query files ###", flush=True)
    # Create the query files
    fasta_count, created_queries = query_files_creation(low_cov_regions, assembly_file, min_query_len,
                                                                       query_files_dir_path, queries_per_file)

    run_information.append("Created query .fasta files:\t\t" + str(fasta_count))
    run_information.append("Created queries:\t\t\t" + str(created_queries))

    # Calculate the runtime
    run_time = time.strftime("%Hh%Mm%Ss", time.gmtime((time.time() - start_time)))
    # Append the runtime info at the second position
    run_information = [run_information[0], str("Runtime:\t\t\t\t" + run_time)] + run_information[1:]

    if verbose_func:
        print("### Create log file ###", flush=True)
    # Create run information file
    run_info_file = open((storage_files_dir_path + "query_creation_" + time.strftime("%Y%m%d-%H%M%S") + ".txt"), "w")
    for line in run_information:
        run_info_file.write((line + "\n"))

    run_info_file.close()

    print("#### clcr.query_creation is succesfully finished! ####", flush=True)

    return None


def main():
    version = "1.0.0"
    # Initialise parser
    parser = argparse.ArgumentParser(description="##### CLCR query creation #####",
                                     epilog="This function detects the regions with low read coverages in a given pbc "
                                            "(per base coverage) file. The sequences of the detected regions are "
                                            "extracted from the handed over assembly and stored as .fasta query files "
                                            "in the query_files directory of the handed over project. The query_files "
                                            "directory will be overwritten!\n"
                                            "A log file with run information and a original_low_cov_regions.tsv with "
                                            "the original detected low coverage regions before the merging step are "
                                            "stored at the storage_files dir.\n"
                                            "ATTENTION!: old query dir and DIAMOND output dir of current project dir is "
                                            "overwritten as preparation for a new cluster run!!")

    # Differentiate between required and optional arguments
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    # Required arguments
    required.add_argument("-p", "--project_dir", action='store', type=str, required=True,
                          help="Path of the project directory")
    required.add_argument("-c", "--cov_file_path", action='store', type=str, required=True,
                          help="Path of the coverage file")
    required.add_argument("-a", "--assembly_file", action='store', type=str, required=True,
                          help="Path of the assembly file")

    # optional arguments
    optional.add_argument("--low_cov_start", action='store', type=int, required=False, default=15,
                          help="Threshold for detecting a low cov region")
    optional.add_argument("--low_cov_end", action='store', type=int, required=False, default=18,
                          help="Threshold for ending a low cov region")
    optional.add_argument("--min_query_len", action='store', type=int, required=False, default=500,
                          help="Minimum query length")
    optional.add_argument("--queries_per_file", action='store', type=int, required=False, default=5000,
                          help="Queries sequences per query file")
    optional.add_argument("--verbose", action='store_true', required=False,
                          help="Run information is print in the command line")

    # Parse args
    args = parser.parse_args()
    # Hand over parsed arguments to create_queries function
    create_queries(args)

