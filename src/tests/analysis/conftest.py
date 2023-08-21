"""
Pytest fixtures for the analysis folder
"""
import pytest

#from tests.conftest import *
from serena.utilities.comparison_structures import (ComparisonNucResults)
from serena.utilities.local_minima_variation import (ComparisonLMV,
                                                    ComparisonLMVResponse)
from serena.analysis.investigator import (SettingsAssertionLMV,
                                          LMVAssertionResult,
                                          RatioResults,
                                          ComparisonEvalResults,
                                          InvestigatorResults)

from serena.analysis.judge_pool import (JudgesResults)



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

@pytest.fixture
def intantiated_investigator_results(initialized_comparison_eval_results:ComparisonEvalResults, comparison_nuc_result:ComparisonNucResults, initialized_comparison_lmv_response:ComparisonLMVResponse, initialized_lmv_assertion_results:LMVAssertionResult):
    return InvestigatorResults(comparison_eval_results=initialized_comparison_eval_results,
                                                                   comp_nuc_counts=comparison_nuc_result,
                                                                   lmv_values=initialized_comparison_lmv_response,
                                                                   lmv_assertions=initialized_lmv_assertion_results,
                                                                   num_groups=3,
                                                                   total_structures_ensemble=5)

"""
Fixtures for judge pool
"""
@pytest.fixture
def initialized_judge_result():
    return JudgesResults(is_good_count=1,
                         is_good_switch=True,
                         is_on_off_count=2,
                         is_on_off_switch=False,
                         is_powerful_count=3,
                         is_powerful_switch=True,
                         switchable_groups_list=[1,2,3],
                         powerfull_groups_list=[2,3,4],
                         on_off_groups_list=[3,4,5])