import pytest
from serena.analysis.investigator import (InvestigatorResults)
from serena.analysis.judge_pool import (JudgesResults,
                                        AnalysisJudgePool)

def test_initialized_judge_result(initialized_judge_result:JudgesResults):
    assert initialized_judge_result.is_good_count == 1
    assert initialized_judge_result.is_good_switch == True
    assert initialized_judge_result.is_on_off_count == 2
    assert initialized_judge_result.is_on_off_switch == False
    assert initialized_judge_result.is_powerful_count == 3
    assert initialized_judge_result.is_powerful_switch == True
    assert initialized_judge_result.switchable_groups_list == [1,2,3]
    assert initialized_judge_result.powerfull_groups_list == [2,3,4]
    assert initialized_judge_result.on_off_groups_list == [3,4,5]

"""
Test judge pool
"""

def test_is_switch_judge(intantiated_investigator_results:InvestigatorResults):
    judge_pool:AnalysisJudgePool = AnalysisJudgePool()
    result:JudgesResults = judge_pool.is_switch_judge(investigator=intantiated_investigator_results)
    assert result.is_good_count==2


