"""
File for class for analysis stuff
"""

from typing import List
from dataclasses import dataclass

from serena.utilities.comparison_structures import ComparisonNucCounts, ComparisonResult
from serena.utilities.ensemble_variation import EV, EVResult
from serena.utilities.local_minima_variation import ComparisonLMV
from serena.utilities.weighted_structures import WeightedNucCounts

@dataclass
class SettingsAssertionLMV():
    diff_limit_mfe:float = 0
    diff_limit_comp:float = 1

@dataclass
class AssertionLMVResult():
    mfe_pronounced:bool = False
    comp_pronounced: bool = False
    rel_pronounced: bool = False

@dataclass
class SwitchabilitySettings():
    limit: float = 1.5 

@dataclass
class SwitchynessResult():
    is_switchable_group:List[bool]
    switchable_groups_list:List[int]
    is_powerfull_switch_group:List[bool]
    powerfull_groups_list:List[int]

@dataclass
class RatioResults():
    unbound_to_total_ratio:float = 0
    bound_ratio: float = 0
    bound_to_both_ratio = 0

@dataclass
class InvestigatorResults():
    ratios: List[RatioResults] 
    comp_nuc_counts: List[ComparisonNucCounts]
    lmv_values: List[ComparisonLMV]
    lmv_assertions: List[AssertionLMVResult]
    num_groups:int = 0

@dataclass
class AnalysisRatioResults():
    last_count_unbound:float=0
    last_count_bound:float=0
    last_count_both: float = 0
    unbound_to_total_ratio:float = 0
    bound_ratio: float = 0
    last_unbound_ratio = 0
    last_bound_ratio = 0
    last_both_ratio = 0
    bound_to_both_ratio = 0

class ComparisonNucAnalysis():

    def __init__(self) -> None:
        pass

    def evaluate_switchability_ensemble(self, comparison_data:List[ComparisonResult], lmv_data:List[ComparisonLMV], lmv_assertions: List[AssertionLMVResult], settings: SwitchabilitySettings):
        
        total_groups:int = len(comparison_data)
        total_nucs: int = comparison_data[0].comp_counts.num_nucs

        is_switchable_group:List[bool] = []
        switchable_groups_list:List[bool] = []
        is_powerfull_switch_group:List[bool] = []
        powerfull_groups_list:List[bool]= []

        is_good_count:int=0
        is_excelent_count:int =0

        is_powerful_switch = False
        is_good_switch = False

        last_unbound:float=0
        last_bound:float=0
        last_both: float = 0

        bound_total_list: List[int] = []
        unbound_total_list: List[int] = []
        

        for group_index in range(total_groups):
            bound: int = comparison_data[group_index].comp_counts.bound_count
            unbound: int= comparison_data[group_index].comp_counts.unbound_count
            both_nuc:int = comparison_data[group_index].comp_counts.bound_count
            dot_nuc:int = comparison_data[group_index].comp_counts.dot_count
            
            unbound_to_total_ratio:float = 0
            bound_ratio: float = 0
            last_unbound_ratio = 0
            last_bound_ratio = 0
            last_both_ratio = 0
            bound_to_both_ratio = 0
            try:
                last_unbound_ratio = last_unbound/unbound 
            except:
                pass
            
            try:
                bound_ratio = bound/unbound
            except:
                pass

            try:

                if bound_hold != -1:
                    #do normal                    
                    if bound_hold < last_bound: 
                        if bound_hold == 0:
                            bound_hold = 1                   
                        last_bound_ratio = bound/bound_hold 
                    else:
                        last_bound_ratio = bound/last_bound 
                else:
                    last_bound_ratio = bound/last_bound

                if bound > last_bound:
                    #its getting bigger so record that
                    bound_hold = last_bound   
                else:
                    bound_hold = -1    
            except:
                pass
            
            try:
                last_both_ratio = both_nuc/last_both 
            except:
                pass
            
            try:
                bound_to_both_ratio = bound/(both_nuc - unbound)
            except:
                pass

            unbound_to_total_ratio = unbound/total_nucs
            bound_to_total_ratio = bound/total_nucs
            both_nuc_total= both_nuc/total_nucs
            dot_nuc_total= dot_nuc/total_nucs

            bound_total_list.append(bound_to_total_ratio)
            unbound_total_list.append(unbound_to_total_ratio)  

            bound_stats: str = f'BURatio:{round(bound_ratio,2)},both_Raise:{round(last_both_ratio,2)} BRaise:{round(last_bound_ratio,2)}, UDrop:{round(last_unbound_ratio,2)},BothTotal:{round(both_nuc_total,2)}, BoundTotal:{round(bound_to_total_ratio,2)}, UTotal:{round(unbound_to_total_ratio,2)}, bound_both:{round(bound_to_both_ratio,2)} B:{bound}, U:{unbound}. both:{both_nuc}'

            limit:float = settings.limit

            last_unbound_ratio = round(last_unbound_ratio,2)
            last_bound_ratio = round(last_bound_ratio,2)
            unbound_to_total_ratio = round(unbound_to_total_ratio,2)
            bound_ratio = round(bound_ratio,2)

            ev_weight_asserted:bool = lmv_assertions[group_index].comp_pronounced
            ev_weigth_under_limit:bool = False
            ev_weight_limit:int = 25
            if lmv_data[group_index].lmv_comp.ev_normalized < ev_weight_limit:
                ev_weigth_under_limit = True 

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
            
            last_unbound = unbound
            last_bound = bound
            last_both = both_nuc

 
            is_switchable_group.append(is_good_switch)
            is_powerfull_switch_group.append(is_powerful_switch)
        
        result: SwitchynessResult = SwitchynessResult(is_switchable_group=is_switchable_group,
                                                      switchable_groups_list=switchable_groups_list,
                                                      is_powerfull_switch_group=is_powerfull_switch_group,
                                                      powerfull_groups_list=powerfull_groups_list)
        return result


class LocalMinimaVariationAnalysis():

    def __init__(self) -> None:
        pass

    def evaluate_lmv_for_structure_presence(self, lmv_results:ComparisonLMV, setting:SettingsAssertionLMV):          

        ev_comp = lmv_results.lmv_comp.ev_normalized
        ev_comp_limit: float = 25
        ev_mfe = lmv_results.lmv_mfe.ev_normalized

        diff_limit_mfe:float = setting.diff_limit_mfe
        diff_limit_comp:float = setting.diff_limit_comp

        lmv_presence_result: AssertionLMVResult = AssertionLMVResult()
               

        diff_comp:float = round(ev_mfe,2) - round(ev_comp,2)
        if round(ev_comp,2) < round(ev_mfe,2) and diff_comp >= diff_limit_comp:
            lmv_presence_result.comp_pronounced = True

        diff_mfe = round(ev_comp,2) - round(ev_mfe,2)
        if round(ev_mfe,2) <= round(ev_comp,2) and (diff_mfe >= diff_limit_mfe):
            lmv_presence_result.mfe_pronounced = True
        
        return lmv_presence_result