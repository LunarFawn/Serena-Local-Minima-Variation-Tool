import pytest
from serena.analysis.investigator import (InvestigatorResults)
from serena.analysis.judge_pool import (JudgesResults,
                                        AnalysisJudgePool,
                                        LMVSwitchJudgeResult,
                                        CompSwitchJudgeResult)

def test_initialized_judge_result(initialized_judge_result:JudgesResults):
    assert initialized_judge_result.comp_switch_judge.is_good_count == 1
    assert initialized_judge_result.comp_switch_judge.is_good_switch == True
    assert initialized_judge_result.comp_switch_judge.is_powerful_count == 3
    assert initialized_judge_result.comp_switch_judge.is_powerful_switch == True
    assert initialized_judge_result.lmv_switch_judge.is_on_off_count == 2
    assert initialized_judge_result.lmv_switch_judge.is_on_off_switch == False
    assert initialized_judge_result.comp_switch_judge.switchable_groups_list == [1,2,3]
    assert initialized_judge_result.comp_switch_judge.powerfull_groups_list == [2,3,4]
    assert initialized_judge_result.lmv_switch_judge.on_off_groups_list == [3,4,5]

"""
Test judge pool
"""

def test_is_comp_switch_judge(intantiated_investigator_results:InvestigatorResults):
    judge_pool:AnalysisJudgePool = AnalysisJudgePool()
    result:CompSwitchJudgeResult = judge_pool.is_comp_switch_judge(investigator=intantiated_investigator_results)
    assert result.is_good_count==0
    assert result.is_good_switch == False
    assert result.is_powerful_count == 0
    assert result.is_powerful_switch == False
    

#def test_is_lmv_switch_judge(intantiated_investigator_results:InvestigatorResults):
#    judge_pool:AnalysisJudgePool = AnalysisJudgePool()
#    result:LMVSwitchJudgeResult = judge_pool.is_lmv_switch_judge(investigator=intantiated_investigator_results)
#    assert result.is_on_off_count == 2
#    assert result.is_on_off_switch == True
#    assert result.on_off_groups_list == [1,2]

#def test_all_judges_function(intantiated_investigator_results:InvestigatorResults):
#    judge_pool:AnalysisJudgePool = AnalysisJudgePool()
#    result: JudgesResults = judge_pool.run_all_judges(investigator=intantiated_investigator_results)
#    assert result.comp_switch_judge.is_good_count==0
#    assert result.comp_switch_judge.is_good_switch == False
#    assert result.comp_switch_judge.is_powerful_count == 0
#    assert result.comp_switch_judge.is_powerful_switch == False
#    assert result.lmv_switch_judge.is_on_off_count == 2
#    assert result.lmv_switch_judge.is_on_off_switch == True
#    assert result.lmv_switch_judge.on_off_groups_list == [1,2]
    