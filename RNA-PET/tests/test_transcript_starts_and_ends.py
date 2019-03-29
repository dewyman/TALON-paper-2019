import pytest
import sys
sys.path.append("..")
import create_intervals as cI
import get_transcript_start_end_intervals as trints

class TestTranscriptInterval(object):
    def test_plus_strand(self):
        """ Interval on the plus strand"""

        dist = 2
        bed_entry = ["chr1", "101", "108", "test_name", "0", "+"]
        start, end = trints.make_intervals(bed_entry, dist) 
     
        assert start == (99, 104)
        assert end == (105, 110)

    def test_minus_strand(self):
        """ Interval on the minus strand """
        
        dist = 2
        bed_entry = ["chr1", "101", "108", "test_name", "0", "-"]
        start, end = trints.make_intervals(bed_entry, dist)

        assert start == (105, 110)
        assert end == (99, 104)
