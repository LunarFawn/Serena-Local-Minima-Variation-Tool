"""
File for class for analysis stuff
"""

from typing import List
from dataclasses import dataclass

from serena.utilities.comparison_structures import ComparisonNucCounts, ComparisonResult, ComparisonNucResults
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
    unbound_to_total_ratio:float
    bound_ratio: float
    last_unbound_ratio: float
    last_bound_ratio: float
    last_both_ratio: float
    bound_to_both_ratio: float
    bound_to_total_ratio:float
    both_nuc_total:float

@dataclass
class ComparisonEvalResults():
    #last_count_unbound:float=0
    #last_count_bound:float=0
    #last_count_both: float = 0
    ratios:List[RatioResults]
    BRaise_list:List[float]
    BUratio_list:List[float]
    bound_total_list: List[int]
    unbound_total_list: List[int]
    nuc_penatly_count:int
    first_BUratio:float

#the code that makes this is not written yet...dont forget
@dataclass
class InvestigatorResults():
    comparison_eval_results: ComparisonEvalResults
    comp_nuc_counts: ComparisonNucResults
    lmv_values: ComparisonLMVResponse
    lmv_assertions: LMVAssertionResult
    num_groups:int = 0
    total_structures_ensemble:int = 0


class ComparisonInvestigator():

    def __init__(self) -> None:
        pass

    def evalulate_comparison_nucs(self, comparison_nucss:ComparisonNucResults)->ComparisonEvalResults:

        BRaise_list:List[float] = []
        BUratio_list:List[float] = []
        bound_total_list: List[int] = []
        unbound_total_list: List[int] = []
        ratios:List[RatioResults] = []
        nuc_penatly_count:int = 0
        for group_index in range(len(comparison_nucss)):
            last_index:int = 0
            if group_index > 0:
                last_index = group_index -1
            unbound:float = comparison_nucss.comparison_nuc_counts[group_index].unbound_count
            last_unbound:float = comparison_nucss.comparison_nuc_counts[last_index].unbound_count
            
            bound:float = comparison_nucss.comparison_nuc_counts[group_index].bound_count
            last_bound:float = comparison_nucss.comparison_nuc_counts[last_index].bound_count

            both_nuc:float = comparison_nucss.comparison_nuc_counts[group_index].both_count
            last_both:float = comparison_nucss.comparison_nuc_counts[last_index].both_count

            dot_nuc:float = comparison_nucss.comparison_nuc_counts[group_index].dot_count

            nuc_count:int = comparison_nucss.comparison_nuc_counts[last_index].num_nucs

            unbound_to_total_ratio:float = 0
            bound_to_total_ratio:float = 0
            both_nuc_total:float = 0
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
            
            #added to address the ones with 0 in the first group
            if group_index > 0:
                if BRaise_list[group_index-1] == 0 and bound > 0:
                    last_bound_ratio = bound


            try:
                last_both_ratio = both_nuc/last_both 
            except:
                pass
            
            try:
                bound_to_both_ratio = bound/(both_nuc - unbound)
            except:
                pass



            unbound_to_total_ratio = unbound/nuc_count
            bound_to_total_ratio = bound/nuc_count
            both_nuc_total= both_nuc/nuc_count
            dot_nuc_total= dot_nuc/nuc_count

            bound_total_list.append(bound_to_total_ratio)
            unbound_total_list.append(unbound_to_total_ratio)  
            
            #now round teh data to make it more managable
            last_unbound_ratio = round(last_unbound_ratio,2)
            last_bound_ratio = round(last_bound_ratio,2)
            unbound_to_total_ratio = round(unbound_to_total_ratio,2)
            bound_ratio = round(bound_ratio,2)
            bound_stats: str = f'BURatio:{round(bound_ratio,2)},both_Raise:{round(last_both_ratio,2)} BRaise:{round(last_bound_ratio,2)}, UDrop:{round(last_unbound_ratio,2)},BothTotal:{round(both_nuc_total,2)}, BoundTotal:{round(bound_to_total_ratio,2)}, UTotal:{round(unbound_to_total_ratio,2)}, bound_both:{round(bound_to_both_ratio,2)} B:{bound}, U:{unbound}. both:{both_nuc}'


            #this is only for the fist kcal group
            if group_index == 0:
                nuc_penatly_count = bound
                first_BUratio = round(bound_ratio,2)
            
            BUratio_list.append(round(bound_ratio,2))
            BRaise_list.append(round(bound,2))

            ratio_results:RatioResults = RatioResults(unbound_to_total_ratio=unbound_to_total_ratio,
                                                      bound_ratio=bound_ratio,
                                                      last_unbound_ratio=last_unbound_ratio,
                                                      last_bound_ratio=last_bound_ratio,
                                                      last_both_ratio=last_both_ratio,
                                                      bound_to_both_ratio=bound_to_both_ratio,
                                                      bound_to_total_ratio=bound_to_total_ratio,
                                                      both_nuc_total=both_nuc_total
                                                      )
            ratios.append(ratio_results)

        comparison_eval_results: ComparisonEvalResults = ComparisonEvalResults(ratios=ratios,
                                                                               BRaise_list=BRaise_list,
                                                                               BUratio_list=BUratio_list,
                                                                               bound_total_list=bound_total_list,
                                                                               unbound_total_list=unbound_total_list,
                                                                               nuc_penatly_count=nuc_penatly_count,
                                                                               first_BUratio=first_BUratio)
        return comparison_eval_results
    


class LocalMinimaVariationInvestigator():

    def __init__(self) -> None:
        pass

    def evaluate_lmv_for_structure_presence(self, lmv_data:ComparisonLMVResponse, setting:SettingsAssertionLMV):          

        ev_comp_limit: float = 25

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