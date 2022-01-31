#!python3

"""This file contains the maintenance functions"""
__author__ = "6947325: Johannes Zieres"
__credits__ = ""
__email__ = "johannes.zieres@gmail.com"


def check_correct_healing(healing_data_tsv_path, healed_assembly_file_path):
    """
    Works only when the scaffolds and the positions in scaffold are sorted ascending
    :param healing_data_tsv_path:
    :param healed_assembly_file_path:
    :return:
    """

    # Read in the assembly file for fast access
    healed_assembly_file = open(healed_assembly_file_path)
    scaffold_list = []  # Filled with the sequences of the fna file
    temp_scaffold = []  # Contains the current scaffold
    current_header = ""  # Containing the header of the current scaffold

    print("Read in assembly file")
    # Filling the .fna region list/ reading in the original assembly
    for line in healed_assembly_file:

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
    healed_assembly_file.close()
    # Delete the first empty initialising element
    scaffold_list.pop(0)

    print("Read in healing information")

    # Read in the healing information, each scaffold separately
    healing_data = []   # Contains sublists for each scaffold, which then contains one sublist for each corrected pos.
    healing_data_file = open(healing_data_tsv_path)
    first_line_split = next(healing_data_file).strip().split()
    current_scaffold_name = first_line_split[0]     # Initialise with first scaffold
    current_scaffold_data = [first_line_split]      # Stores the data of each scaffold, initialised with first line

    # Read in the healing data file
    for line in healing_data_file:
        line_split = line.strip().split()

        if line_split[0] == current_scaffold_name:
            current_scaffold_data.append(line_split)

        # Case for new scaffold
        else:
            healing_data.append(current_scaffold_data)
            current_scaffold_name = line_split[0]
            current_scaffold_data = [line_split]

    # Append last remaining scaffold
    healing_data.append(current_scaffold_data)

    healing_data_file.close()

    correct_position_count = 0      # Counts how many of the N's are correctly placed

    print("Check healing positions")

    temp_false = 0

    for scaffold_healing_positions in healing_data:
        pos_shift_count = 0     # Counts how many Bp the original positions are shifted by N insertions

        current_scaffold_sequence = str     # Initialising
        # Search the relating scaffold
        for scaffold_sequence in scaffold_list:
            if scaffold_sequence[0] == scaffold_healing_positions[0][0]:
                current_scaffold_sequence = scaffold_sequence
                break

        # Sort the healing positions in each scaffold ascending
        sorted_healing_pos = sorted(scaffold_healing_positions, key=lambda temp_healing_pos: int(temp_healing_pos[1]))

        # Checks all healing positions in the current scaffold
        for healing_pos in sorted_healing_pos:

            """ if ([healing_pos[0], healing_pos[3], healing_pos[4]] == ['HiC_scaffold_1', '881155', '883251']) or \
                    ([healing_pos[0], healing_pos[3], healing_pos[4]] == ['HiC_scaffold_1', '1218609', '1219109']):
                interesting_position = int(healing_pos[1]) + pos_shift_count
                print("#################")
                print(healing_pos)
                print("pos_shift_count", pos_shift_count)
                # Gibt das N in der Mitte aus
                print(current_scaffold_sequence[1][(interesting_position - 2):(interesting_position + 3)])
                print(current_scaffold_sequence[1][(interesting_position - 20):(interesting_position + 21)])"""

            if current_scaffold_sequence[1][(int(healing_pos[1]) + pos_shift_count)] == "N" and (temp_false == 0):
                correct_position_count += 1

            # Increase position count
            if healing_pos[2] == "D":
                pos_shift_count += 1
            else:
                pos_shift_count += 2

    print("Correct placed positions: ", correct_position_count)

    return None


def custom_query_files_creation(list_with_scaffold_specific_low_cov_reg_lists, fna_file_path, min_length, query_dir,
                                seq_per_fasta):

    """
    CUSTOM version where the header for each query is saved at the third positions after the start and the end position
    of the query.

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
                    current_fasta.write((">" + scaffold[1][0][2] + "\n"))

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
                    current_fasta.write((">" + scaffold[1][0][2] + "\n"))

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
                    current_fasta.write((">" + scaffold[1][0][2] + "\n"))
                    # Writing the query sequence in the fasta file
                    current_fasta.write(("".join(current_scaffold[1][region_start:region_end]) + "\n"))
                    current_tuple += 1
                    count_added_queries += 1
                    already_added_queries += 1

                break       # Break because there is only one matching scaffold sequence for each scaffold

    # Close the last .fasta file
    current_fasta.close()

    # User information
    print("Amount of created fasta files: ", fasta_count)
    print("Amount of created queries: ", count_added_queries)

    return None


def main():
    print("Maintenance main executed")

    healing_data_tsv_path = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_10.01.2022/" \
                            "storage_files/healing_data.tsv"
    healed_assembly_file_path = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_27.01.2022/" \
                                "healed_assembly/healed_assembly.fna"

    check_correct_healing(healing_data_tsv_path, healed_assembly_file_path)


if __name__ == '__main__':
    main()
