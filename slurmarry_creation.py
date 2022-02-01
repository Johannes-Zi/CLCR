#!python3

"""This file contains the functions for the protein database, slurmjob and queries creation"""
__author__ = "6947325: Johannes Zieres"
__credits__ = ""
__email__ = "johannes.zieres@gmail.com"

import os


def create_slurmarry(protein_database, input_dir, output_dir, slurm_dir):
    """
    Simply creates the slurmarray file for a given input and output directory (each dir with an / at the end!)
    :param protein_database: trivial
    :param input_dir: trivial
    :param output_dir: trivial
    :param slurm_dir: path to the directory where the slurm file should be created
    :return: CLCR_slurmarray.slurm file at the cwd
    """

    # creating the .slurm file for the jobarray
    slurm_file = open(slurm_dir, "w")

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
    # os.system(("sbatch --array=1-" + str(fasta_count) + slurm_filename))

    return None


def main():
    print("database comparison main executed!")


if __name__ == '__main__':
    main()
