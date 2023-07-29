"""
File to hold the code to judge if RNA is a switch
"""

from serena.utilities.comparison_structures import ComparisonNucCounts, ComparisonResult
from serena.utilities.ensemble_variation import EV, EVResult
from serena.utilities.local_minima_variation import ComparisonLMV

class AnalysisJudgePool():

    def __init__(self) -> None:
        pass

    if lmv_data[group_index].lmv_comp.ev_normalized < ev_weight_limit:
            ev_weigth_under_limit = True 

    def functional_judge_001(self):
        
        if (last_unbound_ratio >= limit or last_bound_ratio >= limit) and unbound_to_total_ratio <=.3 and ev_weigth_under_limit is True and bound > 2:
                is_good_switch = True
                switchable_groups_list.append(group_index)
                is_good_count = is_good_count+1


       

            
            
            if last_unbound_ratio >= limit and last_bound_ratio >= limit and bound_ratio >=2 and ev_weight_asserted is True:
                is_powerful_switch = True
                powerfull_groups_list.append(group_index)
                is_excelent_count = is_excelent_count +1

            if (last_unbound_ratio >= limit or last_bound_ratio >= limit) and unbound_to_total_ratio <=.2 and ev_weight_asserted is True:
                is_powerful_switch = True
                powerfull_groups_list.append(group_index)
                is_excelent_count = is_excelent_count +1

            if bound_ratio >=  limit and unbound_to_total_ratio <=.15 and ev_weight_asserted is True:
                is_powerful_switch = True
                powerfull_groups_list.append(group_index)
                is_excelent_count = is_excelent_count +1

            if last_bound_ratio >=  2 and unbound_to_total_ratio <=.2:
                is_powerful_switch = True
                powerfull_groups_list.append(group_index)
                is_excelent_count = is_excelent_count +1
            
            if last_bound_ratio > 3 and ev_weight_asserted is True:
                is_good_switch = True
                switchable_groups_list.append(group_index)
                is_good_count = is_good_count + 1
                is_powerful_switch = True
                powerfull_groups_list.append(group_index)
                is_excelent_count = is_excelent_count + 1