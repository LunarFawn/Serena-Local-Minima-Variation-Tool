import pytest

from serena.analysis.ensemble_analysis import (ReferenceStructures,
                                               ProcessEnsemble,
                                               InvestigateEnsembleResults,
                                               InvestigateEnsemble)
from serena.utilities.ensemble_groups import SingleEnsembleGroup, MultipleEnsembleGroups
from serena.utilities.weighted_structures import WeightedEnsembleResult, WeightedStructure
from serena.utilities.local_minima_variation import ComparisonLMV, ComparisonLMVResponse, LocalMinimaVariation
from serena.utilities.comparison_structures import ComparisonNucCounts, ComparisonResult, ComparisonNucResults

def test_reference_structures(initialized_reference_structure:ReferenceStructures):
    assert initialized_reference_structure.mfe_structure.structure == '(...))'
    assert initialized_reference_structure.weighted_structures.structs[0].structure == '..().)'

def test_process_ensemble_for_weighted_structures(process_ensemble_weighted_result:WeightedEnsembleResult):
    assert process_ensemble_weighted_result.structs[0].structure == '...).)'
    assert process_ensemble_weighted_result.structs[1].structure == '(...))'
    
def test_process_ensemble_for_lmv(multiple_ensemble_groups:MultipleEnsembleGroups, initialized_reference_structure:ReferenceStructures):
    process_ensemble:ProcessEnsemble = ProcessEnsemble()
    result: ComparisonLMVResponse = process_ensemble.process_ensemble_for_lmv(ensemble=multiple_ensemble_groups,
                                                                              ref_structures=initialized_reference_structure)
    assert result.lmv_comps[0].lmv_comp.ev_normalized == 2.0
    assert result.lmv_comps[0].lmv_mfe.ev_normalized == 2.0
    assert result.lmv_comps[0].lmv_rel.ev_normalized == 2.0
    assert result.lmv_comps[1].lmv_comp.ev_normalized == 0.5
    assert result.lmv_comps[1].lmv_mfe.ev_normalized == 4.5
    assert result.lmv_comps[1].lmv_rel.ev_normalized == 0.5

def test_process_ensemble_for_comparison_structures(multiple_ensemble_groups:MultipleEnsembleGroups, process_ensemble_weighted_result:WeightedEnsembleResult):
    process_ensemble:ProcessEnsemble = ProcessEnsemble()
    result: ComparisonNucResults = process_ensemble.process_ensemble_for_comparison_structures(raw_ensemble=multiple_ensemble_groups,
                                                                                               weighted_ensemble=process_ensemble_weighted_result)
    assert result.comparison_nuc_counts[0].both_count == 2
    assert result.comparison_nuc_counts[0].bound_count == 1
    assert result.comparison_nuc_counts[0].dot_count == 0
    assert result.comparison_nuc_counts[0].num_nucs == 6
    assert result.comparison_nuc_counts[0].unbound_count == 3

    assert result.comparison_nuc_counts[1].both_count == 2
    assert result.comparison_nuc_counts[1].bound_count == 4
    assert result.comparison_nuc_counts[1].dot_count == 0
    assert result.comparison_nuc_counts[1].num_nucs == 6
    assert result.comparison_nuc_counts[1].unbound_count == 0

def test_initialized_investigate_ensemble_results(initialized_investigate_ensemble_results:InvestigateEnsembleResults):
    assert initialized_investigate_ensemble_results.advanced_scores.total_score == 6
    assert initialized_investigate_ensemble_results.basic_scores.total_score == 1
    assert initialized_investigate_ensemble_results.number_structures == 10

def test_investigate_and_score_ensemble_simple(multiple_ensemble_groups:MultipleEnsembleGroups):
    investiagtion:InvestigateEnsemble = InvestigateEnsemble()
    result:InvestigateEnsembleResults= investiagtion.investigate_and_score_ensemble(ensemble=multiple_ensemble_groups)
    assert result.basic_scores.bonuses == 4
    assert result.basic_scores.functional_switch_score == 3
    assert result.basic_scores.on_off_switch_score == 1
    assert result.basic_scores.penalties == 0
    assert result.basic_scores.powerful_switch_score == 3
    assert result.basic_scores.total_score == 7
    assert result.advanced_scores.comp_bonus == 0.5
    assert result.advanced_scores.comp_penalty == 0
    assert result.advanced_scores.excess_struct_penalty == 0
    assert result.advanced_scores.lmv_bonus == 0
    assert result.advanced_scores.lmv_penalty == 0
    assert result.advanced_scores.total_score == 0.5
    
    
    