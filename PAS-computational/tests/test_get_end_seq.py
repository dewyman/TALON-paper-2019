import pytest
import sys
import pyfasta as pf
import os
sys.path.append("..")
script_dir = os.path.dirname(os.path.realpath(__file__))
utils_dir = "/".join(script_dir.split("/")[0:-2] + ["RNA-PET"])
sys.path.append(script_dir)
sys.path.append(utils_dir)
import create_intervals as cI
import get_last_n_transcript_seq as glt 

class TestGetEndSeq(object):
    def test_plus_strand(self):
        """ Get sequence at the end of plus strand transcript"""

        entry = ['chr1', '827666', '841738', 'ENST00000441765.5', '0', '+']
        chromosome = entry[0]
        strand = entry[-1]
        dist = 10

        genome_file = "chr1.fa"
        genome = pf.Fasta(genome_file)

        # Get interval for end
        interval_start, interval_end = glt.make_end_interval(entry, dist)
        assert interval_start == 841728
        assert interval_end == 841738

        # Get sequence for end
        seq = glt.fetch_sequence(chromosome, interval_start, interval_end, strand, genome)
        assert seq == "ATGAATATCT"

    def test_minus_strand(self):
        """ Get sequence at the end of minus strand transcript"""

        entry = ['chr1', '185216', '195411', 'ENST00000623083.4', '0', '-']
        chromosome = entry[0]
        strand = entry[-1]
        dist = 35

        genome_file = "chr1.fa"
        genome = pf.Fasta(genome_file) 

        # Get interval for end
        interval_start, interval_end = glt.make_end_interval(entry, dist)
        assert interval_start == 185216
        assert interval_end == 185251 

        # Get sequence for end
        seq = glt.fetch_sequence(chromosome, interval_start, interval_end, strand, genome)
        assert seq == "AGGCAGAGGAGGACGAGGACGACTGGGAATCCTAG"
