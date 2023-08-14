import pytest

from serena.utilities.ensemble_structures import (Sara2SecondaryStructure, 
                                        Sara2StructureList, 
                                        )

from serena.utilities.local_minima_variation import ComparisonLMV, ComparisonLMVResponse
from serena.local_minima_variation import  LocalMinimaVariation
from serena.utilities.ensemble_variation import EV, EVResult
from serena.utilities.ensemble_groups import MultipleEnsembleGroups, SingleEnsembleGroup
from serena.utilities.weighted_structures import WeightedEnsembleResult

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

def test_get_multi_group_lmv(multiple_ensemble_groups:MultipleEnsembleGroups,  secondary_structure_5:Sara2SecondaryStructure):
    local_minima_variation:LocalMinimaVariation = LocalMinimaVariation()
    ev_results:EVResult = local_minima_variation.get_multi_group_lmv(ensemble=multiple_ensemble_groups,
                                                                    reference_structure=secondary_structure_5)
    assert ev_results.ev_values[0].ev_normalized == 3.0
    assert ev_results.ev_values[0].ev_structure == 0
    assert ev_results.ev_values[0].ev_ThresholdNorm == 0
    assert ev_results.ev_values[1].ev_normalized == 0.5
    assert ev_results.ev_values[1].ev_structure == 0
    assert ev_results.ev_values[1].ev_ThresholdNorm == 0

def test_get_single_group_lmv(single_ensemble_group:SingleEnsembleGroup, secondary_structure_5:Sara2SecondaryStructure):
    local_minima_variation:LocalMinimaVariation = LocalMinimaVariation()
    ev_results:EVResult = local_minima_variation.get_single_group_lmv(ensemble_group=single_ensemble_group,
                                                                    reference_structure=secondary_structure_5)
    assert ev_results.ev_values[0].ev_normalized == 3.0
    assert ev_results.ev_values[0].ev_structure == 0
    assert ev_results.ev_values[0].ev_ThresholdNorm == 0

def test_get_relative_multi_group_lmv(multiple_ensemble_groups:MultipleEnsembleGroups):
    local_minima_variation:LocalMinimaVariation = LocalMinimaVariation()
    ev_results:EVResult = local_minima_variation.get_relative_mutli_group_lmv(ensemble=multiple_ensemble_groups)
    assert ev_results.ev_values[0].ev_normalized == 2.0
    assert ev_results.ev_values[0].ev_structure == 0
    assert ev_results.ev_values[0].ev_ThresholdNorm == 0
    assert ev_results.ev_values[1].ev_normalized == 0.5
    assert ev_results.ev_values[1].ev_structure == 0
    assert ev_results.ev_values[1].ev_ThresholdNorm == 0

def test_get_weighted_structure_restul_lmv(multiple_ensemble_groups:MultipleEnsembleGroups, weighted_ensemble_result:WeightedEnsembleResult):
    local_minima_variation:LocalMinimaVariation = LocalMinimaVariation()
    ev_results:EVResult = local_minima_variation.get_weighted_multi_group_lmv(ensemble=multiple_ensemble_groups,
                                                                              weighted_structures=weighted_ensemble_result)
    assert ev_results.ev_values[0].ev_normalized == 2.0
    assert ev_results.ev_values[0].ev_structure == 0
    assert ev_results.ev_values[0].ev_ThresholdNorm == 0
    assert ev_results.ev_values[1].ev_normalized == 0.5
    assert ev_results.ev_values[1].ev_structure == 0
    assert ev_results.ev_values[1].ev_ThresholdNorm == 0