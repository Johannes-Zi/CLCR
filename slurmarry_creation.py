#!python3

"""This file contains the functions for the protein database, slurmjob and queries creation"""
__author__ = "6947325: Johannes Zieres"
__credits__ = ""
__email__ = "johannes.zieres@gmail.com"

import os
import glob
import argparse


def create_slurmarry(protein_database, input_dir, output_dir, slurm_dir):
    """
    Simply creates the slurmarray file for a given input and output directory (each dir with an / at the end!)
    :param protein_database: trivial
    :param input_dir: trivial
    :param output_dir: trivial
    :param slurm_dir: path to the directory where the slurm file should be created
    :return: CLCR_slurmarray.slurm file at the cwd
    """

    slurm_file = slurm_dir + "CLCR_slurmarray.slurm"

    # creating the .slurm file for the jobarray
    slurm_file = open(slurm_file, "w")

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
                     "#SBATCH --job-name=\"CLCR_run\"\n"
                     "#SBATCH --output=\"" + slurm_dir + "slurm-%j.out\"\n\n")

    slurm_file.write(blast_command)
    slurm_file.close()

    # sending the job via slurm to the cluster
    # os.system("sbatch --array=1-" + str(fasta_count) + slurm_filename)

    return None


def prepare_slurm_run(args):
    """
    Creates the slurmarray file for the current run, if auto_run ist set True, the job will be automatically submitted
    to the computer cluster
    :param args: Arguments from parser

    The following parameters are handed over in the form of args
    -> project_dir: path to the project directory
    -> protein_database: path to the protein database
    -> auto_run: automatic job submitting if True

    :return: None
    """

    # Initialise args as values for a better overview
    project_dir = args.project_dir
    protein_database = args.protein_database
    auto_run = args.auto_run

    # Delete the slurm directory with all files in it, if it already exists
    slurm_dir_path = project_dir + "slurm_files/"
    os_command = " if [ -d " + slurm_dir_path + " ] ; then rm -r " + slurm_dir_path + "; fi"
    os.system(os_command)

    # Create the new slurm directory
    os_command = "mkdir " + slurm_dir_path
    os.system(os_command)

    # Relevant dirs for diamond
    input_dir = project_dir + "query_files/"
    output_dir = project_dir + "diamond_output/"

    # Create slurmarray file
    create_slurmarry(protein_database, input_dir, output_dir, slurm_dir_path)

    # Start slurm job array
    if auto_run:
        # Count the input files
        input_file_list = glob.glob(input_dir + "/temp_in_*.fasta")
        file_count = len(input_file_list)

        # Submit job
        os_command = "sbatch --array=1-" + str(file_count) + " " + slurm_dir_path + "CLCR_slurmarray.slurm"
        os.system(os_command)

    return None


def string_to_bool(v):
    """Enables usage of args True/False"""
    if isinstance(v, bool):
        return v
    if v.lower() in ("True", "true"):
        return True
    elif v.lower() in ("False", "false"):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def main():
    version = "1.0.0"
    # Initialise parser
    parser = argparse.ArgumentParser(description="CLCR cluster run",
                                     epilog="The assembly healing with CLCR consists of 3 different steps:\n"
                                            "1. Query creation with the create_queries function\n"
                                            "2. Performing Diamond alignments, when using slurm as scheduler the "
                                            "prepare_slurm_run function could be used\n"
                                            "3. Creating a healed assembly version with the create_healed_assembly "
                                            "function")

    # Differentiate between required and optional arguments
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    # Required arguments
    required.add_argument("-p", "--project_dir", action='store', metavar="", type=str, required=True,
                          help="Path of the project directory")
    required.add_argument("-c", "--protein_database", action='store', metavar="", type=str, required=True,
                          help="Path of the protein database")

    # optional arguments
    optional.add_argument("--auto_run", type=string_to_bool, nargs='?', const=True, default=False,
                          help="Activate automatic slurm job submission.")

    # Parse args
    args = parser.parse_args()
    # Hand over parsed arguments to create_queries function
    prepare_slurm_run(args)

