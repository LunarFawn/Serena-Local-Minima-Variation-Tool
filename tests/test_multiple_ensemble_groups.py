import pytest
from typing import List, Dict, NamedTuple

from serena.utilities.ensemble_structures import (Sara2SecondaryStructure, 
                                        Sara2StructureList, 
                                        KcalRanges)
from serena.utilities.ensemble_groups import SingleEnsembleGroup, MultipleEnsembleGroups
from test_sara_secondary_structure_lists import test_default_new_secondary_struct_list
from test_sara_secondary_structure import test_empty_secondary_struct

def test_empty_multiple_ensemble_groups(empty_multiple_ensemble_groups:MultipleEnsembleGroups):
    assert empty_multiple_ensemble_groups.groups == []
    assert empty_multiple_ensemble_groups.raw_groups == []
    assert empty_multiple_ensemble_groups.non_switch_state_mfe_kcal == 0
    test_empty_secondary_struct(empty_multiple_ensemble_groups.non_switch_state_structure)
    assert empty_multiple_ensemble_groups.switched_state_mfe_kcal == 0
    test_empty_secondary_struct(empty_multiple_ensemble_groups.switched_state_structure)
    assert empty_multiple_ensemble_groups.groups_dict == {}
    assert empty_multiple_ensemble_groups.group_values == []
    assert empty_multiple_ensemble_groups.num_groups == 0
    assert empty_multiple_ensemble_groups.group_kcal_ranges == []

def test_initialized_multiple_ensemble_groups(multiple_ensemble_groups:MultipleEnsembleGroups):
    assert multiple_ensemble_groups.non_switch_state_mfe_kcal == 10
    assert multiple_ensemble_groups.switched_state_mfe_kcal == 20
    assert multiple_ensemble_groups.non_switch_state_structure.structure == '..().)'
    assert multiple_ensemble_groups.switched_state_structure.structure == '(...))'

def test_set_num_groups(empty_multiple_ensemble_groups:MultipleEnsembleGroups):
    empty_multiple_ensemble_groups.num_groups = 10
    assert empty_multiple_ensemble_groups.num_groups == 10

def test_set_groups(empty_multiple_ensemble_groups:MultipleEnsembleGroups ):
    """
    need to make a groups list with a new singel ensemble group
    """
