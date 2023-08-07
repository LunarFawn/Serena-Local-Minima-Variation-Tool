from re import S
import pytest
from typing import List, Dict, NamedTuple

from serena.utilities.ensemble_structures import (Sara2SecondaryStructure, 
                                        Sara2StructureList, 
                                        KcalRanges)
from serena.utilities.ensemble_groups import SingleEnsembleGroup, MultipleEnsembleGroups
from test_sara_secondary_structure_lists import test_default_new_secondary_struct_list
from test_sara_secondary_structure import test_empty_secondary_struct

def test_empty_single_ensemble_group(empty_single_ensemble_group:SingleEnsembleGroup):
    #test that the ensemble list is empty default
    test_default_new_secondary_struct_list(empty_secondary_structure_list=empty_single_ensemble_group.group)
    assert empty_single_ensemble_group.multi_state_mfe_kcal == []
    assert empty_single_ensemble_group.multi_state_mfe_struct == []
    assert empty_single_ensemble_group.kcal_end == 0
    assert empty_single_ensemble_group.kcal_span == 0
    assert empty_single_ensemble_group.kcal_start == 0

def test_set_single_ensemble_group_properties(empty_single_ensemble_group:SingleEnsembleGroup, secondary_structures_list_2_item:Sara2StructureList):
    empty_single_ensemble_group.group = secondary_structures_list_2_item
    assert empty_single_ensemble_group.group.sara_stuctures[0].structure == '((.)))'
    assert empty_single_ensemble_group.group.sara_stuctures[1].structure == '..().)'
    
    mfe_structs_list:List[str] = ['((..))','(...))']
    empty_single_ensemble_group.multi_state_mfe_struct = mfe_structs_list
    assert empty_single_ensemble_group.multi_state_mfe_struct[0] == '((..))'
    assert empty_single_ensemble_group.multi_state_mfe_struct[1] == '(...))'
    
    mfe_kcal_list:List[float] = [-10,-20]
    empty_single_ensemble_group.multi_state_mfe_kcal = mfe_kcal_list
    assert empty_single_ensemble_group.multi_state_mfe_kcal[0] == -10
    assert empty_single_ensemble_group.multi_state_mfe_kcal[1] == -20

    empty_single_ensemble_group.kcal_end = 10
    empty_single_ensemble_group.kcal_span = 20
    empty_single_ensemble_group.kcal_start = 30
    assert empty_single_ensemble_group.kcal_end == 10
    assert empty_single_ensemble_group.kcal_span == 20
    assert empty_single_ensemble_group.kcal_start == 30

def test_fancy_single_ensemble_group_properties(empty_single_ensemble_group:SingleEnsembleGroup):
    empty_single_ensemble_group.append_multi_state_mfe_data('((..))',-10)
    empty_single_ensemble_group.append_multi_state_mfe_data('(...))', -20)
    assert empty_single_ensemble_group.multi_state_mfe_struct[0] == '((..))'
    assert empty_single_ensemble_group.multi_state_mfe_struct[1] == '(...))'
    assert empty_single_ensemble_group.multi_state_mfe_kcal[0] == -10
    assert empty_single_ensemble_group.multi_state_mfe_kcal[1] == -20
    
    empty_single_ensemble_group.update_kcals(start=30,
                                            stop=10,
                                            span=20)
    assert empty_single_ensemble_group.kcal_end == 10
    assert empty_single_ensemble_group.kcal_span == 20
    assert empty_single_ensemble_group.kcal_start == 30

"""
Now multiple ensemble groups
"""

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

def test_set_multiple_ensemble_groups(empty_multiple_ensemble_groups:MultipleEnsembleGroups):
    pass

