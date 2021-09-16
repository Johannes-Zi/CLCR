"""This program is the bruteforce version for detecting interesting regions"""
__author__ = "6947325: Johannes Zieres"
__credits__ = ""
__email__ = "johannes.zieres@gmail.com"


def detect_regions(cov_file_path, cov_start, cov_end):
    """
    Detecting the low coverage reions
    :param cov_file_path: path of input coverage file
    :param cov_start: threshold for detecting a low coverage area
    :param cov_end: threshold for leaving a low coverage area
    :return: list with the low coverage regions
    """

    import_file = open(cov_file_path)  # opens the coverage file
    current_position = 0  # saves the current position
    low_cov = False  # saves the current coverage situation
    low_cov_regions = []  # low coverage regions List with tuple like (start low cov. area, end of low cov. area)
    region_start = int  # saves the start point for a low coverage region

    # detects the low coverage regions and saves them in low_coverage_regions list
    for line in import_file:
        temp_split1 = line.split()

        if low_cov:  # Case for being in a low coverage region
            if int(temp_split1[2]) > cov_end:
                low_cov = False
                low_cov_regions.append((region_start, (current_position - 1)))

        else:  # Case for being in a good coverage region
            if int(temp_split1[2]) <= cov_start:
                low_cov = True
                region_start = current_position

        current_position += 1

    import_file.close()  # closes the input file

    return low_cov_regions


def merge_close_reg(input_list, merge_distance):

    """
    Merging close regions together
    :param input_list: input list, which is containing the low coverage regions
    :param merge_distance: distance threshold, for merging regions which are closer together than the threshold
    :return: region list with the merged regions
    """

    output_list = []
    current_start = int

    current_start = input_list[0][0]    # initialising value

    # merging regions
    for x in range(len(input_list)-1):
        if (input_list[x][1] + merge_distance) < input_list[x+1][0]:
            output_list.append((current_start, input_list[x][1]))
            current_start = input_list[x+1][0]

    # append last region
    # case if the next to last must be merged with the last region and appended
    if (input_list[(len(input_list)-2)][1] + merge_distance) >= input_list[len(input_list)-1][0]:
        output_list.append((current_start, input_list[len(input_list)-1][1]))
    # case for appending a unattached last region
    else:
        output_list.append((input_list[len(input_list)-1][0], input_list[len(input_list)-1][1]))

    return output_list


def combine_regions(input_list_short, input_list_long):

    """
    Searching for combined low coverage regions
    :param input_list_short: list containing the short read low coverage regions
    :param input_list_long: list containing the long read low coverage regions
    :return: list which is containing the regions where both coverages are low
    """

    combined_list = []      # for saving the created combined regions

    # saving the current tuple positions in the lists
    short_list_pos = 0
    long_list_pos = 0

    combined_region = False    # saves the new "combined regions", which are not already represented in the input lists

    # set starting point
    if input_list_short[0][0] <= input_list_long[0][0]:
        startpoint = input_list_short[0][0]
        active_region = "S"     # for short
    else:
        startpoint = input_list_long[0][0]
        active_region = "L"  # for long

    # goes gradually through the regions, for each active region (short or long) the while-loop consideres if the region
    # could be expanded by a region from the other input_list. If a region is expanded the "combined_region" is set TRUE
    # and if the regions could not be extended more the starting point of the region and the current end is appended in
    # combined_list.
    while True:

        # case if the current "active/considered" region is represented by an short read region
        if active_region == "S":

            # break condition
            if (long_list_pos >= len(input_list_long)) or (short_list_pos >= len(input_list_short)):
                break

            # tries every potential region which could expand the current active region
            # while-loop represents every long read region which protrudes in the current active region
            while input_list_long[long_list_pos][0] <= input_list_short[short_list_pos][1]:

                # There are two different possible cases

                # 1. the considered region expands the current active region
                if input_list_long[long_list_pos][1] > input_list_short[short_list_pos][1]:
                    active_region = "L"
                    short_list_pos += 1
                    combined_region = True      # combine argument is only set TRUE if the region is indeed expanded
                    break

                # 2. the considered region is completely in the current active region and could be skipped
                else:
                    long_list_pos += 1
                    if long_list_pos >= len(input_list_long):
                        break

            # case if no region for expanding was found and the current active region represents an end
            if active_region == "S":

                if combined_region:     # the current region is a merged region and will be appended
                    combined_list.append((startpoint, input_list_short[short_list_pos][1]))
                    combined_region = False

                short_list_pos += 1
                if short_list_pos >= len(input_list_short):
                    break

                # setting a new starting point
                if input_list_short[short_list_pos][0] <= input_list_long[long_list_pos][0]:
                    startpoint = input_list_short[short_list_pos][0]
                    active_region = "S"  # for short
                else:
                    startpoint = input_list_long[long_list_pos][0]
                    active_region = "L"  # for long

        # case if the current "active/considered" region is represented by an long read region
        else:

            # break condition
            if (long_list_pos >= len(input_list_long)) or (short_list_pos >= len(input_list_short)):
                break

            # tries every potential region which could expand the current active region
            # while-loop represents every long read region which protrudes in the current active region
            while input_list_short[short_list_pos][0] <= input_list_long[long_list_pos][1]:

                # There are two different possible cases

                # 1. the considered region expands the current active region
                if input_list_short[short_list_pos][1] > input_list_long[long_list_pos][1]:
                    active_region = "S"
                    long_list_pos += 1
                    combined_region = True  # combine argument is only set TRUE if the region is indeed expanded
                    break

                # 2. the considered region is completely in the current active region and could be skipped
                else:
                    short_list_pos += 1
                    if short_list_pos >= len(input_list_short):
                        break

            # case if no region for expanding was found and the current active region represents an end
            if active_region == "L":

                if combined_region:     # the current region is a merged region and will be appended
                    combined_list.append((startpoint, input_list_long[long_list_pos][1]))
                    combined_region = False

                long_list_pos += 1
                if long_list_pos >= len(input_list_long):
                    break

                # setting a new starting point
                if input_list_short[short_list_pos][0] <= input_list_long[long_list_pos][0]:
                    startpoint = input_list_short[short_list_pos][0]
                    active_region = "S"  # for short
                else:
                    startpoint = input_list_long[long_list_pos][0]
                    active_region = "L"  # for long

    # append the last  combined region if the last combined region is at the end of one of the input lists
    if combined_region:
        if active_region == "L":
            combined_list.append((startpoint, input_list_long[long_list_pos][1]))
        else:
            combined_list.append((startpoint, input_list_short[short_list_pos][1]))

    return combined_list


def main():
    print("bruteforce main executed")


if __name__ == '__main__':
    main()
