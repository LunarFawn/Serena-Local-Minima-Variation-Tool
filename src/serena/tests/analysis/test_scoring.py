import pytest

from serena.analysis.scoring import (BasicScoreResults,
                                     AdvancedScoreResults,
                                     SerenaScoring,
                                     JudgesResults) 
from serena.analysis.investigator import (SettingsAssertionLMV,
                                          LMVAssertionResult,
                                          RatioResults,
                                          ComparisonEvalResults,
                                          InvestigatorResults
                                          )

def test_special_penalties_excessive_structs_is_excessive():
    serena_scoring:SerenaScoring = SerenaScoring()
    penaltie: float = serena_scoring.excessive_structures_penalties(10000,2000,5000)
    assert penaltie == 2.5

def test_special_penalties_excessive_structs_inoy_excessive():
    serena_scoring:SerenaScoring = SerenaScoring()
    penaltie: float = serena_scoring.excessive_structures_penalties(10000,2000,10000)
    assert penaltie == 0

def test_empty_basic_score_results(empty_basic_score_results:BasicScoreResults):
    assert empty_basic_score_results.total_score == 0
    assert empty_basic_score_results.functional_switch_score == 0
    assert empty_basic_score_results.powerful_switch_score == 0
    assert empty_basic_score_results.on_off_switch_score == 0
    assert empty_basic_score_results.bonuses == 0
    assert empty_basic_score_results.penalties == 0

def test_initialized_basic_score_results(initialized_basic_score_results:BasicScoreResults):
    assert initialized_basic_score_results.total_score==1
    assert initialized_basic_score_results.functional_switch_score==2
    assert initialized_basic_score_results.powerful_switch_score==3
    assert initialized_basic_score_results.on_off_switch_score==4
    assert initialized_basic_score_results.bonuses==5
    assert initialized_basic_score_results.penalties==6

def test_set_basic_score_results(empty_basic_score_results:BasicScoreResults):
    empty_basic_score_results.bonuses = 1
    empty_basic_score_results.penalties = 2
    empty_basic_score_results.total_score = 3
    empty_basic_score_results.functional_switch_score = 4
    empty_basic_score_results.on_off_switch_score = 5
    empty_basic_score_results.powerful_switch_score = 6
    assert empty_basic_score_results.bonuses == 1
    assert empty_basic_score_results.penalties == 2
    assert empty_basic_score_results.total_score == 3
    assert empty_basic_score_results.functional_switch_score == 4
    assert empty_basic_score_results.on_off_switch_score == 5
    assert empty_basic_score_results.powerful_switch_score == 6

def test_empty_advanced_score_results(empty_advanced_score_results:AdvancedScoreResults):
    assert empty_advanced_score_results.lmv_bonus ==0
    assert empty_advanced_score_results.lmv_penalty == 0    
    assert empty_advanced_score_results.comp_bonus == 0
    assert empty_advanced_score_results.comp_penalty == 0
    assert empty_advanced_score_results.excess_struct_penalty ==0
    assert empty_advanced_score_results.total_score == 0 
    
def test_initialized_advanced_score_results(initialized_advanced_score_results:AdvancedScoreResults):
    assert initialized_advanced_score_results.lmv_bonus == 1
    assert initialized_advanced_score_results.lmv_penalty == 2
    assert initialized_advanced_score_results.comp_bonus == 3
    assert initialized_advanced_score_results.comp_penalty == 4
    assert initialized_advanced_score_results.excess_struct_penalty == 5
    assert initialized_advanced_score_results.total_score == 6

def test_set_advanced_score_results(empty_advanced_score_results:AdvancedScoreResults):
    empty_advanced_score_results.comp_bonus = 1
    empty_advanced_score_results.comp_penalty = 2
    empty_advanced_score_results.excess_struct_penalty = 3
    empty_advanced_score_results.lmv_bonus = 4
    empty_advanced_score_results.lmv_penalty = 5
    empty_advanced_score_results.total_score = 6
    assert empty_advanced_score_results.comp_bonus == 1
    assert empty_advanced_score_results.comp_penalty == 2
    assert empty_advanced_score_results.excess_struct_penalty == 3
    assert empty_advanced_score_results.lmv_bonus == 4
    assert empty_advanced_score_results.lmv_penalty == 5
    assert empty_advanced_score_results.total_score == 6

def test_serena_scoring_basic_score_groups(initialized_judge_result:JudgesResults,intantiated_investigator_results:InvestigatorResults):
    serena_scoring:SerenaScoring = SerenaScoring()
    basic_scores:BasicScoreResults = serena_scoring.basic_score_groups(judge_results=initialized_judge_result,
                                                                       investigator=intantiated_investigator_results)
    assert basic_scores.bonuses == 0
    assert basic_scores.functional_switch_score == 6
    assert basic_scores.on_off_switch_score == 0
    assert basic_scores.penalties == 0
    assert basic_scores.powerful_switch_score == 3
    assert basic_scores.total_score == 9

def test_advanced_score_groups(initialized_judge_result:JudgesResults,intantiated_investigator_results:InvestigatorResults):
    serena_scoring:SerenaScoring = SerenaScoring()
    basic_scores:AdvancedScoreResults = serena_scoring.advanced_score_groups(judge_results=initialized_judge_result,
                                                                          investigator=intantiated_investigator_results)
    assert basic_scores.comp_bonus == 0
    assert basic_scores.comp_penalty == 4
    assert basic_scores.excess_struct_penalty == 0
    assert basic_scores.lmv_bonus == 0
    assert basic_scores.lmv_penalty == 0
    assert basic_scores.total_score == -4
    