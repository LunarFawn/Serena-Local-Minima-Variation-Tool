
import pytest
from typing import List

from serena.utilities.ensemble_structures import Sara2StructureList, Sara2SecondaryStructure
from serena.utilities.weighted_structures import WeightedStructure, WeightedEnsembleResult, WeightedNucCounts, WeightedComparisonResult



def test_weighted_ensemble_result(secondary_structure_3:Sara2SecondaryStructure, secondary_structure_4:Sara2SecondaryStructure):
    struct_list:List[Sara2SecondaryStructure] = []
    struct_list.append(secondary_structure_3)
    struct_list.append(secondary_structure_4)
    result:WeightedEnsembleResult = WeightedEnsembleResult(structs=struct_list)
    assert result.structs[0] == secondary_structure_3
    assert result.structs[1] == secondary_structure_4

def test_weighted_nuc_counts(empty_weighted_nuc_count:WeightedNucCounts):
    assert empty_weighted_nuc_count.num_both == -1
    assert empty_weighted_nuc_count.num_bound == -1
    assert empty_weighted_nuc_count.num_dot == -1
    assert empty_weighted_nuc_count.num_nucs == -1
    assert empty_weighted_nuc_count.num_unbound == -1

def test_initialized_nuc_counts(weighted_nuc_count:WeightedNucCounts):
    assert weighted_nuc_count.num_unbound == 1
    assert weighted_nuc_count.num_both == 2
    assert weighted_nuc_count.num_bound == 3
    assert weighted_nuc_count.num_dot == 4
    assert weighted_nuc_count.num_nucs == 5

def test_setting_weighted_nuc_counts(empty_weighted_nuc_count:WeightedNucCounts):
    empty_weighted_nuc_count.num_both = 1
    empty_weighted_nuc_count.num_bound = 2
    empty_weighted_nuc_count.num_dot = 3
    empty_weighted_nuc_count.num_nucs = 4
    empty_weighted_nuc_count.num_unbound = 5
    assert empty_weighted_nuc_count.num_both == 1
    assert empty_weighted_nuc_count.num_bound == 2
    assert empty_weighted_nuc_count.num_dot == 3
    assert empty_weighted_nuc_count.num_nucs == 4
    assert empty_weighted_nuc_count.num_unbound == 5

def test_empty_weighted_result(empty_weighted_comparison_result:WeightedComparisonResult):
    pass