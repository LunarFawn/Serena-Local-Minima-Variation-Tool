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
                                            InvestigatorResults
                                            )
from serena.analysis.judge_pool import (JudgesResults,
                                        CompSwitchJudgeResult,
                                        LMVSwitchJudgeResult)
from serena.analysis.scoring import (BasicScoreResults,
                                     AdvancedScoreResults,
                                     SerenaScoring)
from serena.utilities.ensemble_structures import (Sara2SecondaryStructure, 
                                        Sara2StructureList, 
                                        KcalRanges)
from serena.utilities.weighted_structures import WeightedEnsembleResult
from serena.analysis.ensemble_analysis import (ReferenceStructures,
                                               ProcessEnsemble,
                                               InvestigateEnsembleResults,
                                               InvestigateEnsemble)
from serena.utilities.ensemble_groups import SingleEnsembleGroup, MultipleEnsembleGroups



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
def initialized_lmv_assertion_results_2():
    return LMVAssertionResult(bound_pronounced=[True, True, False],
                              unbouund_pronounced=[False, False, True],
                              comp_compare_to_mfe=['<', '>', '='],
                              is_on_off_switch=[False,True,True])

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
def realish_ratio_results_1():
    return RatioResults(unbound_to_total_ratio=0.4,
                        both_nuc_total=0.6,
                        bound_ratio=0.5,
                        bound_to_both_ratio=1.0,
                        bound_to_total_ratio=0.2,
                        last_both_ratio=1.0,
                        last_bound_ratio=1.5,
                        last_unbound_ratio=1.0,
                        dot_to_total_ratio=0.8)

@pytest.fixture
def realish_ratio_result_2():
    return RatioResults(unbound_to_total_ratio=0.6,
                        both_nuc_total=0.8,
                        bound_ratio=0.67,
                        bound_to_both_ratio=2.0,
                        bound_to_total_ratio=0.4,
                        last_both_ratio=1.3333333333333333,
                        last_bound_ratio=1.9, #need better test coverage here
                        last_unbound_ratio=0.67,
                        dot_to_total_ratio=1.0)

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
def realish_comparison_eval_results(initialized_ratio_results:RatioResults, realish_ratio_results_1:RatioResults, realish_ratio_result_2:RatioResults):
    return ComparisonEvalResults(ratios=[initialized_ratio_results,realish_ratio_results_1,realish_ratio_result_2],
                                 BRaise_list=[1.1,2.2,3.3],
                                 BUratio_list=[1.2, 2.4,3.5],
                                 bound_total_list=[1,2,5],
                                 unbound_total_list=[2,6,8],
                                 nuc_penatly_count=3,
                                 first_BUratio=5)

@pytest.fixture
def initialized_comparison_eval_results_2(initialized_ratio_results:RatioResults):
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
def intantiated_investigator_results(realish_comparison_eval_results:ComparisonEvalResults, comparison_nuc_result_2:ComparisonNucResults, comparison_lmv_response_yes_switch:ComparisonLMVResponse, initialized_lmv_assertion_results_2:LMVAssertionResult):
    return InvestigatorResults(comparison_eval_results=realish_comparison_eval_results,
                                comp_nuc_counts=comparison_nuc_result_2,
                                lmv_values=comparison_lmv_response_yes_switch,
                                lmv_assertions=initialized_lmv_assertion_results_2,
                                num_groups=3,
                                total_structures_ensemble=5)

"""
Fixtures for judge pool
"""

@pytest.fixture
def comp_swtich_judge_result():
    return CompSwitchJudgeResult(is_good_count=1,
                                    is_good_switch=True,
                                    is_powerful_count=3,
                                    is_powerful_switch=True,
                                    switchable_groups_list=[1,2,3],
                                    powerfull_groups_list=[2,3,4])

@pytest.fixture
def lmv_switch_judge_result():
    return LMVSwitchJudgeResult(is_on_off_count=2,
                                is_on_off_switch=False,
                                on_off_groups_list=[3,4,5])

@pytest.fixture
def initialized_judge_result(comp_swtich_judge_result:CompSwitchJudgeResult, lmv_switch_judge_result:LMVSwitchJudgeResult):
    return JudgesResults(comp_switch_judge=comp_swtich_judge_result,
                         lmv_switch_judge=lmv_switch_judge_result)


"""
fixtures for scoring
"""

@pytest.fixture
def empty_basic_score_results():
    return BasicScoreResults()

@pytest.fixture
def initialized_basic_score_results():
    return BasicScoreResults(total_score=1,
                             functional_switch_score=2,
                             powerful_switch_score=3,
                             on_off_switch_score=4,
                             bonuses=5,
                             penalties=6)

@pytest.fixture
def empty_advanced_score_results():
    return AdvancedScoreResults()

@pytest.fixture
def initialized_advanced_score_results():
    return AdvancedScoreResults(lmv_bonus=1,
                                lmv_penalty=2,
                                comp_bonus=3,
                                comp_penalty=4,
                                excess_struct_penalty=5,
                                total_score=6)

"""
Doing ensemble analysis now
"""

@pytest.fixture
def initialized_reference_structure(secondary_structure_5:Sara2SecondaryStructure, weighted_ensemble_result:WeightedEnsembleResult):
    return ReferenceStructures(mfe_structure=secondary_structure_5,
                               weighted_structures=weighted_ensemble_result)

@pytest.fixture
def process_ensemble_weighted_result(multiple_ensemble_groups:MultipleEnsembleGroups):
    process_ensemble:ProcessEnsemble = ProcessEnsemble()
    result: WeightedEnsembleResult = process_ensemble.process_ensemble_for_weighted_structures(ensemble=multiple_ensemble_groups)
    return result

@pytest.fixture
def initialized_investigate_ensemble_results(initialized_basic_score_results:BasicScoreResults, initialized_advanced_score_results:AdvancedScoreResults):
    return InvestigateEnsembleResults(basic_scores=initialized_basic_score_results,
                                      advanced_scores=initialized_advanced_score_results,
                                      number_structures=10)