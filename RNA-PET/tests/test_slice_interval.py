import pytest
import sys
sys.path.append("..")
import create_intervals as cI

class TestSliceInterval(object):
    def test_dist_zero(self):
        """ Report just the first position on its own (BED), and the last 
            position on its own as zero-based coordinates"""

        transcript_start = 100
        transcript_end = 110
        dist = 0

        # Get interval for start
        interval_start, interval_end = cI.create_interval(transcript_start, 
                                                          "left", dist)
        assert interval_start == 100
        assert interval_end == 101

        # Get interval for end
        interval_start, interval_end = cI.create_interval(transcript_end, 
                                                          "right", dist) 
        assert interval_start == 109
        assert interval_end == 110

    def test_dist_nonzero(self):
        """ Report nonzero interval around the positions """

        transcript_start = 101
        transcript_end = 108
        dist = 2

        # Get interval for start
        interval_start, interval_end = cI.create_interval(transcript_start,
                                                          "left", dist)
        assert interval_start == 99
        assert interval_end == 104

        # Get interval for end
        interval_start, interval_end = cI.create_interval(transcript_end,
                                                          "right", dist)
        assert interval_start == 105
        assert interval_end == 110

class TestCreateEndPieces(object):
    def test_dist_plus_strand(self):
        """ Test plus-strand behavior """
        pos = 110
        strand = "+"
        dist = 2
        interval_start, interval_end = cI.create_end_piece(pos, strand, dist)

        assert interval_start == 108
        assert interval_end == 110

    def test_dist_minus_strand(self):
        """ Test minus-strand behavior """
        pos = 100
        strand = "-"
        dist = 2
        interval_start, interval_end = cI.create_end_piece(pos, strand, dist)

        assert interval_start == 100
        assert interval_end == 102
