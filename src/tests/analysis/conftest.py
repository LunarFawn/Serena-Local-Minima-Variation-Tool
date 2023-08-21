"""
Pytest fixtures for the analysis folder
"""
import pytest

#from tests.conftest import *
from serena.analysis.investigator import (SettingsAssertionLMV,
                                          LMVAssertionResult,
                                          RatioResults,
                                          ComparisonEvalResults,
                                          InvestigatorResults)



@pytest.fixture
def default_settings_assertions_lmv():
    return SettingsAssertionLMV()

@pytest.fixture
def initialized_settings_assertion_lmv():
    return SettingsAssertionLMV(diff_limit_comp=2,
                                diff_limit_mfe=3)
    
@pytest.fixture
def empty_lmv_assertion_result():
    return LMVAssertionResult()

@pytest.fixture
def initialized_lmv_assertion_results():
    return LMVAssertionResult(bound_pronounced=[True, True, False],
                              unbouund_pronounced=[False, False, True],
                              comp_compare_to_mfe=['a', 'b', 'c'],
                              is_on_off_switch=True)

@pytest.fixture
def default_ratio_results():
    return RatioResults()

@pytest.fixture
def initialized_ratio_results():
    return RatioResults(unbound_to_total_ratio=1,
                        both_nuc_total=2,
                        bound_ratio=3,
                        bound_to_both_ratio=4,
                        bound_to_total_ratio=5,
                        last_both_ratio=6,
                        last_bound_ratio=7,
                        last_unbound_ratio=8,
                        dot_to_total_ratio=9)

@pytest.fixture
def empty_comparison_eval_results():
    return ComparisonEvalResults()

@pytest.fixture
def initialized_comparison_eval_results(initialized_ratio_results:RatioResults):
    return ComparisonEvalResults(ratios=[initialized_ratio_results],
                                 BRaise_list=[1.1],
                                 BUratio_list=[1.2],
                                 bound_total_list=[1],
                                 unbound_total_list=[2],
                                 nuc_penatly_count=3,
                                 first_BUratio=5)

@pytest.fixture
def initialized_investigator_results(initialized_comparison_eval_results:ComparisonEvalResults):
    return InvestigatorResults(comparison_eval_results=initialized_comparison_eval_results)