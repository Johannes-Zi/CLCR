#!python3

"""This file contains the functions for the healing of the genome assembly"""
__author__ = "6947325: Johannes Zieres"
__credits__ = ""
__email__ = "johannes.zieres@gmail.com"


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

    print("Read in assembly file")
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

    print("Queries sorted")
    count = 0

    # List with the distance of each healed position to the previous healed position
    healing_position_distribution = []

    for query in sorted_healing_region_list:
        # Search the corresponding scaffold
        print("\033[A                             \033[A")
        print("Current Query", count)

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


def main():
    print("assembly healing main executed!")


if __name__ == '__main__':
    main()