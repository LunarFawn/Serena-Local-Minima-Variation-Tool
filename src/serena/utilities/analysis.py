"""
File for class for analysis stuff
"""

from typing import List
from dataclasses import dataclass

from serena.utilities.comparison_structures import ComparisonNucCounts, ComparisonResult
from serena.utilities.ensemble_variation import EV, EVResult
from serena.utilities.local_minima_variation import ComparisonLMV, ComparisonLMVResponse
from serena.utilities.weighted_structures import WeightedNucCounts

@dataclass
class SettingsAssertionLMV():
    diff_limit_mfe:float = 0
    diff_limit_comp:float = 1

@dataclass
class LMVAssertionResult():
    bound_compare_to_unbound:List[str]
    unbouund_pronounced:List[bool]
    bound_pronounced: List[bool]
    is_on_off_switch:List[bool]

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

#the code that makes this is not written yet...dont forget
@dataclass
class InvestigatorResults():
    ratios: List[RatioResults] 
    comp_nuc_counts: List[ComparisonNucCounts]
    lmv_values: List[ComparisonLMV]
    lmv_assertions: List[LMVAssertionResult]
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

class ComparisonInvestigator():

    def __init__(self) -> None:
        pass
    
    


class LocalMinimaVariationInvestigator():

    def __init__(self) -> None:
        pass

    def evaluate_lmv_for_structure_presence(self, lmv_data:ComparisonLMVResponse, setting:SettingsAssertionLMV):          

        ev_comp = lmv_results.lmv_comp.ev_normalized
        ev_comp_limit: float = 25
        ev_mfe = lmv_results.lmv_mfe.ev_normalized

        diff_limit_mfe:float = setting.diff_limit_mfe
        diff_limit_comp:float = setting.diff_limit_comp

        comp_pronounced:List[bool] = []
        is_on_off_switch:List[bool] = []
        mfe_pronounced:List[bool] = []
        
        for group_index in range(len(lmv_data.lmv_comps)):
            ev_comp:float = lmv_data.lmv_comps[group_index].lmv_comp
            ev_mfe:float = lmv_data.lmv_comps[group_index].lmv_mfe
            
            comp_asserted:bool = False
            is_on_off_:bool = False
            mfe_asserted:bool = False
                

            diff_comp:float = round(ev_mfe,2) - round(ev_comp,2)
            if round(ev_comp,2) < round(ev_mfe,2) and diff_comp >= diff_limit_comp:
                comp_asserted = True
                is_on_off_ = True

            diff_mfe = round(ev_comp,2) - round(ev_mfe,2)
            if round(ev_mfe,2) <= round(ev_comp,2) and (diff_mfe >= diff_limit_mfe):
                mfe_asserted = True
            
            comp_pronounced.append(comp_asserted)
            mfe_pronounced.append(mfe_asserted)
            is_on_off_switch.append(is_on_off_)
        
        ev_comp_to_mfe_list:List[str] = self.bound_comared_unbound_lmv(lmv_data=lmv_data)

        lmv_presence_result: LMVAssertionResult = LMVAssertionResult(bound_compare_to_unbound=ev_comp_to_mfe_list,
                                                                        unbouund_pronounced=mfe_pronounced,
                                                                        bound_pronounced=comp_pronounced,
                                                                        is_on_off_switch=is_on_off_switch)

        return lmv_presence_result
    
    def bound_comared_unbound_lmv(self, lmv_data:ComparisonLMVResponse):
        
        ev_comp_to_mfe_list:List[str] = []

        for group_index in range(len(lmv_data.lmv_comps)):
            ev_comp:float = lmv_data.lmv_comps[group_index].lmv_comp
            ev_mfe:float = lmv_data.lmv_comps[group_index].lmv_mfe
            if ev_comp < ev_mfe:
                ev_comp_to_mfe_list.append('<')
            elif ev_comp == ev_mfe:
                ev_comp_to_mfe_list.append('=')
            elif ev_comp > ev_mfe:
                ev_comp_to_mfe_list.append('>')
        
        return ev_comp_to_mfe_list