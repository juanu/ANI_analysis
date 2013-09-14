#!/usr/local/bin/python
#Created on 4/2/2013
__author__ = 'Juan A. Ugalde'


def split_sequence_fragments(sequence, block_size):
    """
    Function that takes a sequence of letters (usually a DNA sequence), and returns
    fragments of a defined block size
    """
    fragments = []
    total_size = len(sequence)

    for i in range(0, total_size, block_size):
        fragments.append(sequence[i:i + block_size])

    return fragments


def run_blastn(folder, reference, query, names):
    """
    Function that takes a reference file and a query file and run blastn.
    The required input is a folder to create the temporal blast database, and
    reference and query files in fasta format.
    The results will be saved with the name of query versus reference
    """
    import os
    #Make blastdb
    if not os.path.isfile(reference):
        print "Reference file: %s not found" % reference

    #os.system('formatdb -i %s -p F -n %s/reference' % (reference, folder))
    os.system('makeblastdb -in %s -dbtype nucl -out %s/reference' % (reference, folder))

    num_processors = 4  # Number of processors to use for Blast

    query_name, reference_name = names

    blast_output_name = folder + "/" + query_name + "_" + reference_name

    #os.system('blastall -p blastn -a %d -d %s/reference -i %s -X 150 -q -1 -F F -m 8 -o %s' %
    #          (num_processors, folder, query, blast_output_name))

    os.system('blastn -num_threads %d -db %s/reference -query %s -xdrop_gap 150  -penalty -1 -dust no -outfmt 6 '
              '-gapopen 5 -gapextend 2 -out %s' %
    (num_processors, folder, query, blast_output_name))

    return blast_output_name


def get_blast_top_hit(blast_file):
    """
    Parse the blast file. Select the top hit
    """
    blast_results = [line.rstrip() for line in open(blast_file)]
    blast_top_hit = {}

    for blast_line in blast_results:
        best_hit = True
        (queryId, subjectId, percIdentity, alnLength, mismatchCount, gapOpenCount, queryStart,
         queryEnd, subjectStart, subjectEnd, evalue, bitScore) = blast_line.split("\t")

        #get the top hit
        if queryId in blast_top_hit:
            if float(bitScore) < float(blast_top_hit.get(queryId)[11]):
                best_hit = False

        if best_hit:
            blast_top_hit[queryId] = blast_line.split("\t")

    return blast_top_hit


def calculate_ani(blast_results, fragment_length):
    """
    Takes the input of the blast results, and calculates the ANI versus the reference genome
    """
    sum_identity = float(0)
    number_hits = 0  # Number of hits that passed the criteria
    total_aligned_bases = 0  # Total of DNA bases that passed the criteria
    total_unaligned_fragments = 0
    total_unaligned_bases = 0

    conserved_dna_bases = 0

    for query in blast_results:
        identity = blast_results[query][2]
        queryEnd = blast_results[query][7]
        queryStart = blast_results[query][6]

        perc_aln_length = (float(queryEnd) - float(queryStart)) / fragment_length[query]

        if float(identity) > float(69.9999) and float(perc_aln_length) > float(0.69999):
            sum_identity += float(identity)
            number_hits += 1
            total_aligned_bases += fragment_length[query]

        else:
            total_unaligned_fragments += 1
            total_unaligned_bases += fragment_length[query]

        if float(identity) > float(89.999):
            conserved_dna_bases += fragment_length[query]

    return sum_identity, number_hits, total_aligned_bases, total_unaligned_fragments, total_unaligned_bases


def average_ani_results(ani_dictionary):
    """
    This function takes the dictionary that contains the ani dictionary, take the reference and query
    and takes the average between the two results of the combination of reference and query
    """
    refined_ani_results = {}

    for pair in ani_dictionary:
        reference_query_value = ani_dictionary[pair]
        reference, query = pair
        query_reference_value = ani_dictionary[(query, reference)]

        average_value = (reference_query_value + query_reference_value) / 2

        if (query, reference) in refined_ani_results:
            continue
        else:
            refined_ani_results[pair] = average_value

    return refined_ani_results


def create_distance_matrix(ani_dictionary):
    """

    """
    from itertools import count
    import numpy as np

    data = []
    for pair in ani_dictionary:
        reference, query = pair
        value = 100 - float(ani_dictionary[pair])
        data.append([reference, query, value])
        data.append([query, reference, value])

    rows = dict(zip(sorted(set(line[0] for line in data)), count()))
    cols = dict(zip(sorted(set(line[1] for line in data)), count()))

    ani_array = np.zeros((len(rows), len(rows)), dtype=float)

    for row, col, val in data:
        index = (rows[row], cols[col])
        ani_array[index] = val

    return rows, cols, ani_array


if __name__ == '__main__':
    import sys
    import shutil
    import argparse
    import os
    import itertools
    from Bio import SeqIO
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import scipy.cluster.hierarchy as sch
    import scipy.spatial.distance

    program_description = "Script that takes a list of genomes, containing the location of the fasta files," \
                          "and generates a matrix with the ANI values for all the combinations"

    parser = argparse.ArgumentParser(description=program_description)

    parser.add_argument("-i", "--genome_input_list", type=str, help="List with the genome names and files",
                        required=True)

    parser.add_argument("-o", "--output_directory", type=str, help="Output directory", required=True)

    args = parser.parse_args()

    #Create output directory
    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)

    #Create temporal folder for the blast analysis
    temp_folder = args.output_directory + "/temp"
    if not os.path.exists(temp_folder):
        os.makedirs(temp_folder)

    #Read the genome list:
    genome_info = {element[0]: element[1] for element in [line.split("\t") for line in [line.rstrip() for line in open(args.genome_input_list)] if line.strip()]}

    #Create log file
    log_output = open(args.output_directory + "/logfile.txt", 'w')
    mapping_summary = open(args.output_directory + "/mapping_summary.txt", 'w')

    mapping_summary.write("Reference\tReference Genome size\t"
                          "Query\tQuery Genome Size\t"
                          "Total number fragments\tMapped Fragments\tIdentity\t"
                          "Mapped Bases\tUnmapped fragments\tUnmapped Bases\n")

    #Parameters for blast and fragments
    fragment_size = 500

    #Create genome combinations for blast analysis
    genome_combinations = itertools.permutations(genome_info.keys(), 2)
    genome_pair_identity = {}  # Results

    raw_ani_results = {}  # Results of the ANI analysis

    for genome_pair in genome_combinations:
        reference, query = genome_pair[0], genome_pair[1]

        reference_file = genome_info[reference]
        query_file = genome_info[query]

        #Check that the files exists
        if not os.path.isfile(reference_file):
            print "The reference fasta for %s was not found" % reference
            sys.exit("Check the path for the files")

        if not os.path.isfile(query_file):
            print "The query fasta for %s was not found" % query
            sys.exit("Check the path for the files")

        #Create query file, with fragments of 500bp

        query_fragments_file = open(temp_folder + "/query.fna", 'w')

        fragment_number = 1  # Id of each fragment
        genome_query_fragments = 0
        fragment_length_dict = {}  # Store the size of each fragment
        complete_query_genome_size = 0  # Total size of the query genome
        trimmed_query_genome_size = 0  # Total size of genome no Ns

        for seq_record in SeqIO.parse(query_file, "fasta"):
            genome_sequence = seq_record.seq
            edited_genome_sequence = (str(genome_sequence)).replace("N", "")

            fragments = split_sequence_fragments(edited_genome_sequence, fragment_size)
            complete_query_genome_size += len(seq_record.seq)
            trimmed_query_genome_size += len(edited_genome_sequence)

            genome_query_fragments += len(fragments)

            for fragment in fragments:

                fragment_name = "Fragment" + str(fragment_number)
                query_fragments_file.write(">" + fragment_name + "\n" + str(fragment) + "\n")

                fragment_length_dict[fragment_name] = len(fragment)

                fragment_number += 1

        query_fragments_file.close()

        #Print total number of fragments
        log_output.write("For the query genome: %s \n" % query)
        log_output.write("Genome size: %d \n" % complete_query_genome_size)
        log_output.write("Genome size, with no Ns: %d\n" % trimmed_query_genome_size)
        log_output.write("Number of fragments: %d \n" % genome_query_fragments)

        fragment_query_file = temp_folder + "/query.fna"

        #Print information to screen
        sys.stderr.write("Running blast of %s versus %s \n" % (reference, query))
        sys.stderr.flush()

        #Run blast
        blast_file = run_blastn(temp_folder, reference_file, fragment_query_file, ("reference", "query"))

        #Parse the blast result
        blast_top_hit = get_blast_top_hit(blast_file)

        sum_identity, number_hits, total_aligned_bases, total_unaligned_fragments, total_unaligned_bases = \
            calculate_ani(blast_top_hit, fragment_length_dict)

        try:
            reference_query_ani = sum_identity / number_hits
        except ZeroDivisionError:  # Cases were there are no hits
            reference_query_ani = 0

        #Store the results
        raw_ani_results[(reference, query)] = reference_query_ani

        #Get the size of the reference genome
        reference_genome_size = 0
        for seq_record in SeqIO.parse(reference_file, "fasta"):
            reference_genome_size += len(seq_record.seq)

        results = [reference, str(reference_genome_size), query, str(trimmed_query_genome_size),
                   str(genome_query_fragments), str(number_hits), str(reference_query_ani),
                   str(total_aligned_bases), str(total_unaligned_fragments), str(total_unaligned_bases)]

        mapping_summary.write("\t".join(results) + "\n")

    ##Take the average of the reference query values

    final_ani_results = average_ani_results(raw_ani_results)

    #Generate matrix file
    rows, cols, ani_array = create_distance_matrix(final_ani_results)
    order_col_labels = sorted(cols, key=cols.get)

    #Save matrix file
    matrix_file = open(args.output_directory + "/matrix_file.txt", 'w')
    matrix_file.write("\t" + "\t".join(order_col_labels) + "\n")

    for row_label, row in zip(order_col_labels, ani_array):
        matrix_file.write(row_label + "\t" + "\t".join(str(n) for n in row) + "\n")

    #Run hierarchical analysis and save the plot

    distance_matrix = scipy.spatial.distance.squareform(ani_array)
    linkage_matrix = sch.linkage(distance_matrix, method="complete", metric="euclidean")  # Method and metric
    X = sch.dendrogram(linkage_matrix, labels=order_col_labels, orientation="left")

    plt.subplots_adjust(left=0.3)
    plt.savefig(args.output_directory + "/ANI_hier_plot.pdf")

    #Close final files
    log_output.close()
    mapping_summary.close()
    matrix_file.close()

    #Remove the temporal folder
    shutil.rmtree(temp_folder)