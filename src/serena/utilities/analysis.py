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
class LMVAssertionResult():
    bound_compare_to_unbound:List[str]
    unbouund_pronounced:bool
    bound_pronounced: bool
    is_on_off_switch:bool

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

    def evaluate_lmv_for_structure_presence(self, lmv_results:ComparisonLMV, setting:SettingsAssertionLMV):          

        ev_comp = lmv_results.lmv_comp.ev_normalized
        ev_comp_limit: float = 25
        ev_mfe = lmv_results.lmv_mfe.ev_normalized

        diff_limit_mfe:float = setting.diff_limit_mfe
        diff_limit_comp:float = setting.diff_limit_comp

        comp_pronounced:bool = False
        is_on_off_switch:bool = False
        mfe_pronounced:bool = False
               

        diff_comp:float = round(ev_mfe,2) - round(ev_comp,2)
        if round(ev_comp,2) < round(ev_mfe,2) and diff_comp >= diff_limit_comp:
            comp_pronounced = True
            is_on_off_switch = True

        diff_mfe = round(ev_comp,2) - round(ev_mfe,2)
        if round(ev_mfe,2) <= round(ev_comp,2) and (diff_mfe >= diff_limit_mfe):
            mfe_pronounced = True
        

        lmv_presence_result: LMVAssertionResult = LMVAssertionResult(bound_compare_to_unbound=ev_comp_to_mfe,
                                                                     unbouund_pronounced=mfe_pronounced,
                                                                     bound_pronounced=comp_pronounced,
                                                                     is_on_off_switch=is_on_off_switch)

        return lmv_presence_result
    
    def bound_comared_unbound_lmv(self, ev_comp:float, ev_mfe:float):
        
        ev_comp_to_mfe:List[str] = []
        if ev_comp < ev_mfe:
            ev_comp_to_mfe.append('<')
        elif ev_comp == ev_mfe:
            ev_comp_to_mfe.append('=')
        elif ev_comp > ev_mfe:
            ev_comp_to_mfe.append('>')
            