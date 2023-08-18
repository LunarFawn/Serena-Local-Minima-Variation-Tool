"""
File to hold the code to judge if RNA is a switch
"""

from dataclasses import dataclass
from typing import List

from serena.utilities.comparison_structures import ComparisonNucCounts, ComparisonResult
from serena.utilities.ensemble_variation import EV, EVResult
from serena.utilities.local_minima_variation import ComparisonLMV
from src.serena.analysis.investigator import InvestigatorResults

@dataclass
class SwitchabilitySettings():
    """
    Class for holding values for limits when
    deterining switchability
    """
    limit: float = 1.5 

@dataclass
class JudgesResults():
    """
    Class to hold the results from
    the judge decisions
    """
    is_good_switch:bool
    is_powerful_switch:bool
    is_on_off_switch:bool

    is_good_count:int
    is_powerful_count:int
    is_on_off_count: int 
    
    switchable_groups_list:List[int]
    powerfull_groups_list:List[int]
    on_off_groups_list:List[int]

class AnalysisJudgePool():
    """
    Class for all the different specialized judges
    """
    def __init__(self) -> None:
        pass
        #is_powerful_switch:bool = False
        #is_good_switch:bool = False
        #is_good_count:int = 0
        #is_excelent_count:int = 0
        #current_group_index:int = -1

    def is_switch_judge(self, investigator:InvestigatorResults):
                
        num_groups: int = investigator.num_groups
        
        limit: float = 1.5 
        is_switchable_group:List[bool] = []
        switchable_groups_list:List[int] = []
        is_powerfull_switch_group:List[bool] = []
        powerfull_groups_list:List[int] = []
        is_good_count:int = 0
        is_excelent_count:int = 0
        is_powerful_switch:bool = False
        is_good_switch:bool = False

        for current_group_index in range(len(current_group_index)):
            last_index:int = 0
            if current_group_index>1:
                last_index = current_group_index-1

            last_unbound_ratio:float = investigator.ratios[last_index].unbound_to_total_ratio
            last_unbound_ratio = round(last_unbound_ratio,2)
            last_bound_ratio: float = investigator.ratios[last_index].bound_ratio
            last_bound_ratio = round(last_bound_ratio,2)        
            unbound_to_total_ratio:float = investigator.ratios[current_group_index].unbound_to_total_ratio
            unbound_to_total_ratio = round(unbound_to_total_ratio,2)     
            bound_ratio:float = investigator.ratios[current_group_index].bound_ratio
            bound_ratio = round(bound_ratio,2)

            bound: int = investigator.comp_nuc_counts[current_group_index].bound_count

            lmv_data:List[ComparisonLMV] = investigator.lmv_values
            ev_weight_asserted:bool = investigator.lmv_assertions[current_group_index].comp_pronounced
            ev_weigth_under_limit:bool = False
            ev_weight_limit:int = 25
            if lmv_data[current_group_index].lmv_comp.ev_normalized < ev_weight_limit:
                ev_weigth_under_limit = True 

            if (last_unbound_ratio >= limit or last_bound_ratio >= limit) and unbound_to_total_ratio <=.3 and ev_weigth_under_limit is True and bound > 2:
                is_good_switch = True
                #switchable_groups_list.append(current_group_index)
                is_good_count = is_good_count+1

                
            if last_unbound_ratio >= limit and last_bound_ratio >= limit and bound_ratio >=2 and ev_weight_asserted is True:
                is_powerful_switch = True
                powerfull_groups_list.append(current_group_index)
                is_excelent_count = is_excelent_count +1

            if (last_unbound_ratio >= limit or last_bound_ratio >= limit) and unbound_to_total_ratio <=.2 and ev_weight_asserted is True:
                is_powerful_switch = True
                powerfull_groups_list.append(current_group_index)
                is_excelent_count = is_excelent_count +1

            if bound_ratio >=  limit and unbound_to_total_ratio <=.15 and ev_weight_asserted is True:
                is_powerful_switch = True
                powerfull_groups_list.append(current_group_index)
                is_excelent_count = is_excelent_count +1

            if last_bound_ratio >=  2 and unbound_to_total_ratio <=.2:
                is_powerful_switch = True
                powerfull_groups_list.append(current_group_index)
                is_excelent_count = is_excelent_count +1

            if last_bound_ratio > 3 and ev_weight_asserted is True:
                is_good_switch = True
                is_powerful_switch = True
                is_good_count = is_good_count + 1
                is_excelent_count = is_excelent_count + 1
                switchable_groups_list.append(current_group_index)
                powerfull_groups_list.append(current_group_index)


        results: JudgesResults = JudgesResults(is_excelent_count=is_excelent_count,
                                               is_good_count=is_good_count,
                                               is_good_switch=is_good_switch,
                                               is_powerful_switch=is_powerful_switch,
                                               switchable_groups_list=switchable_groups_list,
                                               powerfull_groups_list=powerfull_groups_list)

        return results
