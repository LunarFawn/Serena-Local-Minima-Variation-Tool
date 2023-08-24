import pytest
from typing import List

from serena.utilities.comparison_structures import (ComparisonNucResults)
from serena.utilities.local_minima_variation import (ComparisonLMV,
                                                    ComparisonLMVResponse)
from serena.analysis.investigator import (SettingsAssertionLMV,
                                            LMVAssertionResult,
                                            RatioResults,
                                            ComparisonEvalResults,
                                            InvestigatorResults,
                                            ComparisonInvestigator,
                                            LocalMinimaVariationInvestigator)



def test_default_settings_assertion(default_settings_assertions_lmv:SettingsAssertionLMV):
    """
    Test SettingsAssertionLMV class with default values
    """
    assert default_settings_assertions_lmv.diff_limit_comp == 1
    assert default_settings_assertions_lmv.diff_limit_mfe == 0

def test_initialized_settings_assertion_lmv(initialized_settings_assertion_lmv:SettingsAssertionLMV):
    """
    Test SettingsAssertionLMV class with initializized values
    """
    assert initialized_settings_assertion_lmv.diff_limit_comp == 2
    assert initialized_settings_assertion_lmv.diff_limit_mfe == 3

def test_empty_lmv_assertion_result(empty_lmv_assertion_result:LMVAssertionResult):
    assert empty_lmv_assertion_result.comp_compare_to_mfe == []
    assert empty_lmv_assertion_result.bound_pronounced == []
    assert empty_lmv_assertion_result.is_on_off_switch == []
    assert empty_lmv_assertion_result.unbouund_pronounced == []

def test_initialized_lmv_assertion_results(initialized_lmv_assertion_results:LMVAssertionResult):
    assert initialized_lmv_assertion_results.bound_pronounced == [True, True, False]
    assert initialized_lmv_assertion_results.unbouund_pronounced==[False, False, True]
    assert initialized_lmv_assertion_results.comp_compare_to_mfe == ['a', 'b', 'c']
    assert initialized_lmv_assertion_results.is_on_off_switch==True

def test_default_ratio_results(default_ratio_results:RatioResults):
    assert default_ratio_results.unbound_to_total_ratio == -1
    assert default_ratio_results.bound_ratio == -1
    assert default_ratio_results.last_unbound_ratio == -1
    assert default_ratio_results.last_bound_ratio == -1
    assert default_ratio_results.last_both_ratio == -1
    assert default_ratio_results.bound_to_both_ratio == -1
    assert default_ratio_results.bound_to_total_ratio == -1
    assert default_ratio_results.both_nuc_total == -1
    assert default_ratio_results.dot_to_total_ratio == -1

def test_initialized_ratio_results(initialized_ratio_results:RatioResults):
    assert initialized_ratio_results.unbound_to_total_ratio == 1
    assert initialized_ratio_results.both_nuc_total == 2
    assert initialized_ratio_results.bound_ratio == 3
    assert initialized_ratio_results.bound_to_both_ratio == 4
    assert initialized_ratio_results.bound_to_total_ratio == 5
    assert initialized_ratio_results.last_both_ratio == 6
    assert initialized_ratio_results.last_bound_ratio == 7
    assert initialized_ratio_results.last_unbound_ratio == 8
    assert initialized_ratio_results.dot_to_total_ratio == 9

def test_empty_comparison_eval_results(empty_comparison_eval_results:ComparisonEvalResults):
    assert empty_comparison_eval_results.ratios == []
    assert empty_comparison_eval_results.BRaise_list ==[]
    assert empty_comparison_eval_results.BUratio_list == []
    assert empty_comparison_eval_results.bound_total_list == []
    assert empty_comparison_eval_results.unbound_total_list == []
    assert empty_comparison_eval_results.nuc_penatly_count == 0
    assert empty_comparison_eval_results.first_BUratio == 0

def test_initialized_comparison_eval_results(initialized_comparison_eval_results:ComparisonEvalResults):
    assert initialized_comparison_eval_results.ratios == [RatioResults(unbound_to_total_ratio=1,
                                                                        both_nuc_total=2,
                                                                        bound_ratio=3,
                                                                        bound_to_both_ratio=4,
                                                                        bound_to_total_ratio=5,
                                                                        last_both_ratio=6,
                                                                        last_bound_ratio=7,
                                                                        last_unbound_ratio=8,
                                                                        dot_to_total_ratio=9)]
    assert initialized_comparison_eval_results.BRaise_list == [1.1]
    assert initialized_comparison_eval_results.BUratio_list == [1.2]
    assert initialized_comparison_eval_results.bound_total_list == [1]
    assert initialized_comparison_eval_results.unbound_total_list == [2]
    assert initialized_comparison_eval_results.nuc_penatly_count == 3
    assert initialized_comparison_eval_results.first_BUratio == 5
    
def test_instantiated_investigator_results(initialized_comparison_eval_results:ComparisonEvalResults, comparison_nuc_result_2:ComparisonNucResults, comparison_lmv_response_yes_switch:ComparisonLMVResponse, initialized_lmv_assertion_results_2:LMVAssertionResult):
    investigator_results:InvestigatorResults = InvestigatorResults(comparison_eval_results=initialized_comparison_eval_results,
                                                                   comp_nuc_counts=comparison_nuc_result_2,
                                                                   lmv_values=comparison_lmv_response_yes_switch,
                                                                   lmv_assertions=initialized_lmv_assertion_results_2,
                                                                   num_groups=3,
                                                                   total_structures_ensemble=5)
    assert investigator_results.comparison_eval_results == initialized_comparison_eval_results
    assert investigator_results.comp_nuc_counts == comparison_nuc_result_2
    assert investigator_results.lmv_values == comparison_lmv_response_yes_switch
    assert investigator_results.lmv_assertions == initialized_lmv_assertion_results_2
    assert investigator_results.num_groups == 3
    assert investigator_results.total_structures_ensemble == 5

"""
Test the comparison investigator
"""

def test_evalulate_comparison_nucs(comparison_nuc_result:ComparisonNucResults):
    investigator:ComparisonInvestigator = ComparisonInvestigator()
    result:ComparisonEvalResults = investigator.evalulate_comparison_nucs(comparison_nucs=comparison_nuc_result)
    assert result.ratios == [RatioResults(unbound_to_total_ratio=0.4,
                                            both_nuc_total=0.6,
                                            bound_ratio=0.5,
                                            bound_to_both_ratio=1.0,
                                            bound_to_total_ratio=0.2,
                                            last_both_ratio=1.0,
                                            last_bound_ratio=1.0, #need better test coverage here
                                            last_unbound_ratio=1.0,
                                            dot_to_total_ratio=0.8),
                            RatioResults(unbound_to_total_ratio=0.6,
                                        both_nuc_total=0.8,
                                        bound_ratio=0.67,
                                        bound_to_both_ratio=2.0,
                                        bound_to_total_ratio=0.4,
                                        last_both_ratio=1.3333333333333333,
                                        last_bound_ratio=2.0, #need better test coverage here
                                        last_unbound_ratio=0.67,
                                        dot_to_total_ratio=1.0)] 
    assert result.BRaise_list == [1,2] 
    assert result.BUratio_list == [0.5, 0.67] 
    assert result.bound_total_list == [0.2,0.4] 
    assert result.unbound_total_list == [0.4,0.6] 
    assert result.nuc_penatly_count == 1
    assert result.first_BUratio == 0.5

"""
Test the Local Minima Investigator 
"""

def test_bound_compared_unbound_lmv(comparison_lmv_response_no_switch:ComparisonLMVResponse):
    investigator:LocalMinimaVariationInvestigator = LocalMinimaVariationInvestigator()
    result:List[str] = investigator.comp_compared_mfe_lmv(lmv_data=comparison_lmv_response_no_switch)
    assert result[0] == '<'
    assert result[1] == '>'
    assert result[2] == '='
    
def test_evaluate_lmv_for_structure_presence_no_lmv_switch(comparison_lmv_response_no_switch:ComparisonLMVResponse,default_settings_assertions_lmv:SettingsAssertionLMV):
    investigator:LocalMinimaVariationInvestigator = LocalMinimaVariationInvestigator()
    result:LMVAssertionResult = investigator.evaluate_lmv_for_structure_presence(lmv_data=comparison_lmv_response_no_switch,
                                                                                setting=default_settings_assertions_lmv)
    assert result.is_on_off_switch == [False,False,False]
    assert result.comp_compare_to_mfe == ['<','>','=']
    assert result.bound_pronounced == [True,False,False]
    assert result.unbouund_pronounced == [False, True, True]

def test_evaluate_lmv_for_structure_presence_yes_lmv_switch(comparison_lmv_response_yes_switch:ComparisonLMVResponse,default_settings_assertions_lmv:SettingsAssertionLMV):
    investigator:LocalMinimaVariationInvestigator = LocalMinimaVariationInvestigator()
    result:LMVAssertionResult = investigator.evaluate_lmv_for_structure_presence(lmv_data=comparison_lmv_response_yes_switch,
                                                                                setting=default_settings_assertions_lmv)
    assert result.is_on_off_switch == [False,True,False]
    assert result.comp_compare_to_mfe == ['>','<','=']
    assert result.bound_pronounced == [False,True,False]
    assert result.unbouund_pronounced == [True, False, True]
