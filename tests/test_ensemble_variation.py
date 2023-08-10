import pytest

from serena.utilities.ensemble_variation import EVResult, EV, EnsembleVariation, EV_Token
from serena.utilities.ensemble_structures import Sara2SecondaryStructure, Sara2StructureList


def test_empty_ev(empty_ev:EV):
    assert empty_ev.ev_normalized == -1
    assert empty_ev.ev_structure == -1
    assert empty_ev.ev_ThresholdNorm == -1

def test_initialized_ev(initialized_ev:EV):
    assert initialized_ev.ev_normalized == 1.1
    assert initialized_ev.ev_structure == 2.2
    assert initialized_ev.ev_ThresholdNorm == 3.3

def test_set_ev_values(empty_ev:EV):
    empty_ev.ev_normalized = 2
    empty_ev.ev_structure = 4
    empty_ev.ev_ThresholdNorm = 6
    assert empty_ev.ev_normalized == 2
    assert empty_ev.ev_structure == 4
    assert empty_ev.ev_ThresholdNorm == 6   

def test_ensemble_variation_algorithm(secondary_structures_list_2_item: Sara2StructureList, secondary_structure_5:Sara2SecondaryStructure):
    ensemble_variation:EnsembleVariation = EnsembleVariation()
    result:EV = ensemble_variation.ensemble_variation_algorithm(kcal_group_structures_list=secondary_structures_list_2_item,
                                                    ref_structure=secondary_structure_5)
    assert result.ev_normalized == 3.0
    assert result.ev_structure == 0
    assert result.ev_ThresholdNorm == 0

def test_ev_result(ev_result:EVResult):
    assert ev_result.ev_values[0].ev_normalized == 1.1
    assert ev_result.ev_values[0].ev_structure == 2.2
    assert ev_result.ev_values[0].ev_ThresholdNorm == 3.3
    assert ev_result.ev_values[1].ev_normalized == 4.4
    assert ev_result.ev_values[1].ev_structure == 5.5
    assert ev_result.ev_values[1].ev_ThresholdNorm == 6.6

def test_empty_3_group_ev_token_group_results(empty_ev_token_3_groups:EV_Token):
    #first test group results initialization
    assert len(empty_ev_token_3_groups.group_results) == 3
    assert empty_ev_token_3_groups.group_results[0].ev_normalized == -1
    assert empty_ev_token_3_groups.group_results[0].ev_structure == -1
    assert empty_ev_token_3_groups.group_results[0].ev_ThresholdNorm == -1
    assert empty_ev_token_3_groups.group_results[1].ev_normalized == -1
    assert empty_ev_token_3_groups.group_results[1].ev_structure == -1
    assert empty_ev_token_3_groups.group_results[1].ev_ThresholdNorm == -1
    assert empty_ev_token_3_groups.group_results[2].ev_normalized == -1
    assert empty_ev_token_3_groups.group_results[2].ev_structure == -1
    assert empty_ev_token_3_groups.group_results[2].ev_ThresholdNorm == -1

def test_empty_3_group_ev_token_group_dict(empty_ev_token_3_groups:EV_Token):
    #first test group results initialization
    assert len(empty_ev_token_3_groups.group_dict) == 3
    assert empty_ev_token_3_groups.group_dict[0].ev_normalized == -1
    assert empty_ev_token_3_groups.group_dict[0].ev_structure == -1
    assert empty_ev_token_3_groups.group_dict[0].ev_ThresholdNorm == -1
    assert empty_ev_token_3_groups.group_dict[1].ev_normalized == -1
    assert empty_ev_token_3_groups.group_dict[1].ev_structure == -1
    assert empty_ev_token_3_groups.group_dict[1].ev_ThresholdNorm == -1
    assert empty_ev_token_3_groups.group_dict[2].ev_normalized == -1
    assert empty_ev_token_3_groups.group_dict[2].ev_structure == -1
    assert empty_ev_token_3_groups.group_dict[2].ev_ThresholdNorm == -1

def test_empty_3_group_ev_token_group_values(empty_ev_token_3_groups:EV_Token):
    #first test group results initialization
    assert len(empty_ev_token_3_groups.group_values) == 3
    assert empty_ev_token_3_groups.group_values[0] == ''   
    assert empty_ev_token_3_groups.group_values[1] == ''
    assert empty_ev_token_3_groups.group_values[2] == ''

def test_empty_3_group_ev_token_group_done_status(empty_ev_token_3_groups:EV_Token):
    #first test group results initialization
    assert len(empty_ev_token_3_groups.group_done_status) == 3
    assert empty_ev_token_3_groups.group_done_status[0] == False  
    assert empty_ev_token_3_groups.group_done_status[1] == False
    assert empty_ev_token_3_groups.group_done_status[2] == False
    assert empty_ev_token_3_groups.is_done == False 




