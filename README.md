# CLCR

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## The CLCR approach
To improve the completeness of draft genome assemblies we created the tool CLCR, which stands for
Correction of Low Coverage Regions. The tool provides functions for short indel detection
and correction in high coverage genome assemblies, in the form of error detection in typically 
error prone assembly regions with a significantly lower read coverage. The detection
of frameshifts is implemented by the integration of the alignment tool Diamond in combination with 
the usage of a protein database with the sequences of closely related organisms.
The frameshift detection by Diamond is followed by the localization of the frameshifts in the
underlying genome in combination with extensive filter heuristics to avoid false correction.
Reading frames that are disturbed by frameshift mutations, putatively caused by a sequencing error,
are locally healed by the insertion of N’s in the assembly.

## Installation
CLCR can be retrieved as a PyPI package, it requires a Python version of 3.7 or higher and a 
[Diamond](https://github.com/bbuchfink/diamond) installation is mandatory. 

The program can be installed with the following command: <br />
```
# Install pip, if nescessary
sudo apt update
sudo apt install python3-pip

# Update pip
pip install --upgrade pip

# Install CLCR
pip install clcr
```
If there are problems with that, the code can be directly cloned from github:
```
# Install git, if nescessary
sudo apt install git-all

# Clone CLCR code from github
git clone https://github.com/Johannes-Zi/CLCR.git

# Move to the cloned directory
pip install .
# Or manually add the program to PATH when there are problems with pip
```
## Quickstart
### Query creation with *clcr.query_creation*
The first step is the creation of the query sequences for the Diamond 
blastx runs. For this, the function detects the regions with low read coverages in a given pbc (per base coverage) file. The
sequences of the detected regions are extracted from the handed over assembly and stored as .fasta query
files in the query_files directory of the handed over project. The query_files directory will be
overwritten! A log file with run information and a original_low_cov_regions.tsv with the original detected
low coverage regions before the merging step are stored at the storage_files dir. <br /> 
*ATTENTION!*: old query directory and DIAMOND output directory of current project dir is overwritten as preparation for a new cluster run!
<br /> 
<br /> 
There are three mandatory arguments required:
* genome assembly
* the matching per base coverage file
* project directory for the CLCR run

There are the following options:
```
  -h, --help            show the help message and exit

required arguments:
  -p PROJECT_DIR, --project_dir PROJECT_DIR
                        Path of the project directory
                        
  -c COV_FILE_PATH, --cov_file_path COV_FILE_PATH
                        Path of the coverage file
                        
  -a ASSEMBLY_FILE, --assembly_file ASSEMBLY_FILE
                        Path of the assembly file

optional arguments:
  --low_cov_start LOW_COV_START
                        Threshold for detecting a low cov region
                        
  --low_cov_end LOW_COV_END
                        Threshold for ending a low cov region
                        
  --min_query_len MIN_QUERY_LEN
                        Minimum query length
                        
  --queries_per_file QUERIES_PER_FILE
                        Queries sequences per query file
                        
  --verbose             Run information is print in the command line
```

## Diamond cluster run



## Contact
For questions, comments or suggestions contact us via [email](mailto:johannes.zieres@gmail.com)
