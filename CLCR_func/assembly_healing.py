#!python3

"""This file contains the functions for the healing of the genome assembly"""
__author__ = "6947325: Johannes Zieres"
__credits__ = ""
__email__ = "johannes.zieres@gmail.com"

import time
import os
import argparse
import CLCR_func.output_processing as output_processing
import CLCR_func.query_creation as query_creation


def heal_assembly_file(healing_region_list, input_fna_path, outut_dir, verbose_func):
    """
    Gets the healing_region_list which contains the frameshift positions in each query, and inserts N's at those
    positions, to heal the reading frame. The healed assembly is afterwards saved in the give output directory.
    The new assembly file doesnt contain comments.
    :param healing_region_list: Contains for each query: scaffold, start pos., end pos. in scaff., frameshift pos. list
    :param outut_dir: directory where the new "healed assembly" is stored (should contain a \ at the last position)
    :param input_fna_path:  File path to the original assembly file
    :param verbose_func: verbose bool
    :return: File path to the new modified assembly file
    """

    # Read in the .fna file for fast access
    input_fna_file = open(input_fna_path)
    scaffold_list = []  # Filled with the sequences of the fna file
    temp_scaffold = []  # Contains the current scaffold
    current_header = ""  # Containing the header of the current scaffold

    if verbose_func:
        print("## Read in original assembly ##", flush=True)
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

    # Correct the frameshifts
    # Inserting the N's from the end to the beginning, to prevent the necessity of adapting the insertion position in
    # dependency to the previous inserted frameshifts
    # For this, in this first sort, the queries in each scaffold are sorted descending by their start position
    # but the reverse option also sorts the scaffolds descending
    sorted_healing_region_list = sorted(healing_region_list, reverse=True,  key=lambda temp_query: (temp_query[0],
                                                                                                    int(temp_query[1])))

    count = 0

    # List with the distance of each healed position to the previous healed position
    healing_position_distribution = []

    if verbose_func:
        print("## Create healed assembly version ##", flush=True)
        print("#", len(sorted_healing_region_list), "Queries with healing positions #", flush=True)
        print(" ")      # empty line for first line clear
    for query in sorted_healing_region_list:
        # Search the corresponding scaffold
        if verbose_func:
            print("\033[A                             \033[A", flush=True)
            print("Current Query", count, flush=True)

        # Saves the position in scaffold of the previous healed frameshift, initialised new for each
        previous_healing_position = None

        for scaffold in scaffold_list:
            # Case where matching scaffold is found
            if scaffold[0] == query[0]:
                query_start_pos = int(query[1])
                count += 1
                # Second sort step, now the frameshift positions in each query are sorted descending
                temp_sorted_query = sorted(query[3], reverse=True,  key=lambda temp_query: int(temp_query[0]))

                for frameshift_position in temp_sorted_query:      # Same usage of reversed() as before
                    # Deletion case
                    if frameshift_position[1] == "D":
                        scaffold[1] = scaffold[1][:(query_start_pos + frameshift_position[0])] + ["N"] + \
                                      scaffold[1][(query_start_pos + frameshift_position[0]):]
                    # Insertion case
                    else:
                        scaffold[1] = scaffold[1][:(query_start_pos + frameshift_position[0])] + ["N"] + ["N"] + \
                                      scaffold[1][(query_start_pos + frameshift_position[0]):]

                    # If for the case where the first position in scaffold was already appended (so for all following
                    # positions in the current scaffold)
                    if previous_healing_position:
                        calculated_distance = previous_healing_position - (query_start_pos + frameshift_position[0])
                        healing_position_distribution.append(calculated_distance)
                        previous_healing_position = (query_start_pos + frameshift_position[0])

                    # Initialise value with first position of new scaffold
                    else:
                        previous_healing_position = (query_start_pos + frameshift_position[0])
                    # Calculate the distance of the current healing position to the previous healed position

                break

    simplified_distance_distribution = [[0, 0], [5, 0], [10, 0], [25, 0], [50, 0], [250, 0], [500, 0], [1000, 0],
                                        [5000, 0], [float("inf"), 0]]

    # Assign each distance to a distance "range" for a better overview
    for calculated_distance in healing_position_distribution:

        # Search the matching "range-group"
        for distance_range in simplified_distance_distribution:
            if calculated_distance <= distance_range[0]:
                distance_range[1] += 1
                break

    # Create the new assembly .fna file path, located in the same dir as the original assembly file
    new_fna_file_path = outut_dir + "healed_assembly.fna"

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

    if verbose_func:
        print("## Healed assembly version saved ##", flush=True)

    return new_fna_file_path, simplified_distance_distribution


def create_detailed_healing_information_file(considered_diamond_hits_list, output_file_path):
    """
    This function simply saves the data which is stored in the considered_diamond_hits_list in a .tsv file.
    This means storing the relevant information for each healed position in the assembly.
    Each line represents one corrected frameshift, the columns are structured like this:
    Scaffold, frameshift pos. in scaff., deletion = D/insertion = I, underlying query start pos. in scaff.,
    query end pos. in scaff., underlying low coverage region start pos. in scaff, low cov. end pos. in scaff.,
    protein hit, e-value, bit-score, similarity-percentage

    Note, the healing positions in each scaffold are sorted ascending, but the sorting of the scaffolds depends on how
    the python sort function behaves with the used scaffold names!

    With this information, the underlying diamond hit for each healing positions could be clearly identified in the
    diamond output data, for further analysis
    :param considered_diamond_hits_list: Output of the filter_out_relevant_results function
    :param output_file_path: path of the output .tsv file, containing the file name
    :return:
    """

    # Create the lines for the output file
    output_line_data_list = []      # Stores the data for the lines

    # Iterate through the query information and save the information in the new file
    for query_data in considered_diamond_hits_list:

        sorted_query_data = sorted(query_data[9], key=lambda temp_query_data: int(temp_query_data[0]))

        for healing_pos in sorted_query_data:

            output_line_data_list.append([query_data[0], str(int(healing_pos[0]) + int(query_data[3])),
                                          healing_pos[1], query_data[3], query_data[4], query_data[1],
                                          query_data[2], query_data[5], query_data[6], query_data[7],
                                          query_data[8]])

    # Sort the healing positions in the list ascending by their healing positions in the scaffolds
    output_lines_list = sorted(output_line_data_list, key=lambda temp_query_data: (temp_query_data[0],
                                                                                   int(temp_query_data[1])))

    # Create new output file
    output_tsv_file = open(output_file_path, "w")

    for out_line_data in output_lines_list:
        output_tsv_file.write(("\t".join(out_line_data) + "\n"))

    output_tsv_file.close()

    return None


def create_healed_assembly(args):
    """
    This function creates a healed version of the handed over assembly. Log file is saved in the storage_files directory
    :param args: Arguments from parser

    The following parameters are handed over in the form of args:
    -> project_dir: path of the project dir
    -> unhealed_assembly: path to the original assembly
    -> dynamic_threshold_dist: The max_detect_distance defines the distance from
                            a detected frameshift position to the original low cov. region, where a frameshift is still
                            considered and not excluded in the further analysis. Values around 10 might be reasonable.
    :return: None
    """

    # Initialise args as values for a better overview
    project_dir = args.project_dir
    unhealed_assembly = args.unhealed_assembly
    dynamic_threshold_dist = args.dynamic_threshold_dist
    verbose_func = args.verbose

    if verbose_func:
        print("#### clcr.assembly_healing called! ####", flush=True)

    # Stores the relevant data of the current run
    run_information = ["CLCR create_healed_assembly run \t\t" + time.ctime(time.time())]   # Initialising

    start_time = time.time()

    if verbose_func:
        print("### Read in Diamond output ###", flush=True)
    # Read in the diamond results
    output_dir = project_dir + "diamond_output/"
    all_diamond_results = output_processing.read_in_diamond_output(output_dir)

    # Create storage files dir if not present
    storage_files_dir_path = project_dir + "storage_files/"
    os_command = " if ! [ -d " + storage_files_dir_path + " ] ; then mkdir " + storage_files_dir_path + "; fi"
    os.system(os_command)

    if verbose_func:
        print("### Read in original low cov regions ###", flush=True)
    # Read in the original low cov. regions
    low_cov_storage_tsv = project_dir + "storage_files/original_low_cov_regions.tsv"
    stored_low_cov_regions = query_creation.read_in_low_cov_tsv_file(low_cov_storage_tsv)

    if verbose_func:
        print("### Filter out relevant frameshifts ###", flush=True)
    # Filter out relevant frameshift positions
    temp_list_1 = output_processing.filter_out_relevant_results(all_diamond_results, dynamic_threshold_dist,
                                                                stored_low_cov_regions)
    considered_diamond_hits_list, healing_region_list, considered_frameshifts_count, frameshifts_detected, \
    insertion_count, deletion_count = temp_list_1

    if verbose_func:
        print("### Create healing data file ###", flush=True)
    # Create healing data file
    healing_data_path = project_dir + "storage_files/healing_data.tsv"
    create_detailed_healing_information_file(considered_diamond_hits_list, healing_data_path)

    # Create healed assembly dir if not present
    storage_files_dir_path = project_dir + "healed_assembly/"
    os_command = " if ! [ -d " + storage_files_dir_path + " ] ; then mkdir " + storage_files_dir_path + "; fi"
    os.system(os_command)

    if verbose_func:
        print("### Heal assembly file ###", flush=True)
    # Create healed assembly file
    new_assembly_dir = project_dir + "healed_assembly/"
    temp_list_2 = heal_assembly_file(healing_region_list, unhealed_assembly, new_assembly_dir, verbose_func)
    new_fna_file_path, simplified_distance_distribution = temp_list_2

    # Calculate the runtime
    run_time = time.strftime("%Hh%Mm%Ss", time.gmtime((time.time() - start_time)))
    # Append relevant information
    run_information.append("Runtime:\t\t\t\t\t" + run_time)
    run_information.append("Queries with at least one frameshift: \t\t" + str(len(healing_region_list)))
    run_information.append("Queries with at least one Diamond hit: \t\t" + str(len(all_diamond_results)))
    run_information.append("Considered frameshifts in the healing process: \t" + str(considered_frameshifts_count))
    # Regardless if they are in the original low cov or not, putative intron transition frameshift excluded
    # run_information.append("Considered frameshifts by the overlapping heuristic: \t" + str(frameshifts_detected))
    run_information.append("Healed insertions: \t\t\t\t" + str(insertion_count))
    run_information.append("Healed deletion: \t\t\t\t" + str(deletion_count))

    if verbose_func:
        print("### Create log file ###", flush=True)
    # Create run information file
    storage_files_dir_path = project_dir + "storage_files/"
    run_info_file = open((storage_files_dir_path + "create_healed_assembly" + time.strftime("%Y%m%d-%H%M%S") + ".txt"),
                         "w")
    for line in run_information:
        run_info_file.write((line + "\n"))

    run_info_file.close()

    if verbose_func:
        print("#### clcr.assembly_healing finished! ####", flush=True)

    return None


def main():
    version = "1.0.0"
    # Initialise parser
    parser = argparse.ArgumentParser(description="##### CLCR query creation #####",
                                     epilog="Function for the creation of a healed assembly version. For this the "
                                            "detected frameshifts in the DIAMOND blastx output are evaluated, "
                                            "extensively filtered and used to created a adapted assembly version with "
                                            "locally healed reading frames.\n"
                                            "The healed assembly version is stored in the healed_assembly dir, and log "
                                            "file for the CLCR run is stored in the storage_files dir.")

    # Differentiate between required and optional arguments
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    # Required arguments
    required.add_argument("-p", "--project_dir", action='store', type=str, required=True,
                          help="Path of the project directory")
    required.add_argument("-c", "--unhealed_assembly", action='store', type=str, required=True,
                          help="Path of the original unhealed assembly file")

    # optional arguments
    optional.add_argument("--dynamic_threshold_dist", action='store', type=int, required=False, default=10,
                          help="The max_detect_distance defines the distance from a detected frameshift position to"
                               " the original low cov. region, where a frameshift is still considered and not excluded"
                               " in the further analysis.")
    optional.add_argument("--verbose", action='store_true', required=False,
                          help="Run information is print in the command line")

    # Parse args
    args = parser.parse_args()
    # Hand over parsed arguments to create_queries function
    create_healed_assembly(args)
