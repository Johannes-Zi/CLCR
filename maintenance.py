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


def main():
    print("Maintenance main executed")

    healing_data_tsv_path = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_10.01.2022/" \
                            "storage_files/healing_data.tsv"
    healed_assembly_file_path = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_27.01.2022/" \
                                "healed_assembly/healed_assembly.fna"

    check_correct_healing(healing_data_tsv_path, healed_assembly_file_path)


if __name__ == '__main__':
    main()
