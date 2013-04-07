Author: Juan A. Ugalde
juanuu@gmail.com

First version, April 4th, 2013.

This is a script to obtain the average nucleotide identity between a group of genomes, using the method of
Goris *et al.*, 2007 (http://ijsb.sgmjournals.org/content/57/1/81.short).

In brief, the method consist in using taking a pair of genomes (reference and query), fragment the query genome in
fragments of size 500bp., and use Blastn.

This script requires Python 2.7 or higher and the libraries:
- Scipy www.scipy.org
- Biopython www.biopython.org

The easiest way to use this script (and other Python related software) is to use a Python distribution that
already includes all the required libraries, like Anaconda http://continuum.io/downloads.html

###Using the script

The input for this script is a tabular file where the first column is the name of the genome, and the second column
is the path for the nucleotide sequence (either single or multiple fasta).
The second required parameter is the name of the output folder where the results will be stored

####Output files:
mapping_summary.txt, a summary of the analysis of each genome
logfile.txt, the query genome, size and number of fragments. Useful to keep track of the status of the program

matrix_file.txt, the identity matrix between all the compared genomes. It indicates the degree of dissimilarity between genomes:
100%-ANI value.

ANI_hier_plot.pdf, a hierarchical plot of the matrix file.



