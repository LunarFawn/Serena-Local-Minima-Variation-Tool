import pytest

from serena.utilities.ensemble_structures import (Sara2SecondaryStructure, 
                                        Sara2StructureList, 
                                        )

from serena.utilities.local_minima_variation import ComparisonLMV, ComparisonLMVResponse
from serena.utilities.ensemble_variation import EV
from serena.utilities.ensemble_groups import MultipleEnsembleGroups

def test_empty_comparison_lmv(empty_comparison_lmv:ComparisonLMV):
    assert empty_comparison_lmv.lmv_comp.ev_normalized == -1
    assert empty_comparison_lmv.lmv_comp.ev_structure == -1
    assert empty_comparison_lmv.lmv_comp.ev_ThresholdNorm == -1
    assert empty_comparison_lmv.lmv_mfe.ev_normalized == -1
    assert empty_comparison_lmv.lmv_mfe.ev_structure == -1
    assert empty_comparison_lmv.lmv_mfe.ev_ThresholdNorm == -1
    assert empty_comparison_lmv.lmv_rel.ev_normalized == -1
    assert empty_comparison_lmv.lmv_rel.ev_structure == -1
    assert empty_comparison_lmv.lmv_rel.ev_ThresholdNorm == -1

def test_initiailized_comparison_lmv(initiailized_comparison_lmv:ComparisonLMV):
    assert initiailized_comparison_lmv.lmv_comp.ev_normalized == 1
    assert initiailized_comparison_lmv.lmv_comp.ev_structure == 2
    assert initiailized_comparison_lmv.lmv_comp.ev_ThresholdNorm == 3
    assert initiailized_comparison_lmv.lmv_mfe.ev_normalized == 4
    assert initiailized_comparison_lmv.lmv_mfe.ev_structure == 5
    assert initiailized_comparison_lmv.lmv_mfe.ev_ThresholdNorm == 6
    assert initiailized_comparison_lmv.lmv_rel.ev_normalized == 7
    assert initiailized_comparison_lmv.lmv_rel.ev_structure == 8
    assert initiailized_comparison_lmv.lmv_rel.ev_ThresholdNorm == 9

def test_set_comparison_lmv(empty_comparison_lmv:ComparisonLMV):
    empty_comparison_lmv.lmv_comp=EV(ev_normalized=1,
                                     ev_structure=2,
                                     ev_ThresholdNorm=3)
    empty_comparison_lmv.lmv_mfe=EV(ev_normalized=4,
                                    ev_ThresholdNorm=6,
                                    ev_structure=5)
    empty_comparison_lmv.lmv_rel=EV(ev_normalized=7,
                                   ev_structure=8,
                                   ev_ThresholdNorm=9)
    assert empty_comparison_lmv.lmv_comp.ev_normalized == 1
    assert empty_comparison_lmv.lmv_comp.ev_structure == 2
    assert empty_comparison_lmv.lmv_comp.ev_ThresholdNorm == 3
    assert empty_comparison_lmv.lmv_mfe.ev_normalized == 4
    assert empty_comparison_lmv.lmv_mfe.ev_structure == 5
    assert empty_comparison_lmv.lmv_mfe.ev_ThresholdNorm == 6
    assert empty_comparison_lmv.lmv_rel.ev_normalized == 7
    assert empty_comparison_lmv.lmv_rel.ev_structure == 8
    assert empty_comparison_lmv.lmv_rel.ev_ThresholdNorm == 9

"""
LMV stuff
"""

def test_get_multi_group_lmv(multiple_ensemble_groups:MultipleEnsembleGroups,  secondary_structure_5:Sara2SecondaryStructure)
    pass