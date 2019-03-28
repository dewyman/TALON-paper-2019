import pytest
import sys
sys.path.append("..")
import create_intervals as cI
import get_RNA_PET_starts_and_ends as grpets

class TestSliceInterval(object):
    def test_plus_strand(self):
        """ Interval on the plus strand"""

        bed_entry = ["chr1", "101", "108", "test_name", "0", "+"]
        start, end = grpets.make_intervals(bed_entry) 
     
        assert start == (101, 102)
        assert end == (107, 108)

    def test_minus_strand(self):
        """ Interval on the minus strand """
        
        bed_entry = ["chr1", "101", "108", "test_name", "0", "-"]
        start, end = grpets.make_intervals(bed_entry)

        assert start == (107, 108)
        assert end == (101, 102)
