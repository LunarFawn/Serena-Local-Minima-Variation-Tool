import pytest

from serena.utilities.comparison_structures import (
        ComparisonStructures, 
        ComparisonResult, 
        ComparisonNucCounts, 
        ComparisonNucResults
)

def test_empty_comparison_nuc_count(empty_comparison_nuc_count:ComparisonNucCounts):
    assert empty_comparison_nuc_count.bound_count == -1
    assert empty_comparison_nuc_count.unbound_count == -1
    assert empty_comparison_nuc_count.both_count == -1
    assert empty_comparison_nuc_count.dot_count == -1
    assert empty_comparison_nuc_count.num_nucs == -1

def test_populated_comparison_nuc_count(comparison_nuc_count:ComparisonNucCounts):
    assert comparison_nuc_count.bound_count == 1
    assert comparison_nuc_count.unbound_count == 2
    assert comparison_nuc_count.both_count == 3
    assert comparison_nuc_count.dot_count == 4
    assert comparison_nuc_count.num_nucs == 5

def test_set_populated_comparison_nuc_count(empty_comparison_nuc_count:ComparisonNucCounts):
    empty_comparison_nuc_count.bound_count = 2
    empty_comparison_nuc_count.unbound_count = 4
    empty_comparison_nuc_count.both_count = 6
    empty_comparison_nuc_count.dot_count = 8
    empty_comparison_nuc_count.num_nucs  = 10
    assert empty_comparison_nuc_count.bound_count == 2
    assert empty_comparison_nuc_count.unbound_count == 4
    assert empty_comparison_nuc_count.both_count == 6
    assert empty_comparison_nuc_count.dot_count == 8
    assert empty_comparison_nuc_count.num_nucs == 10