"""
File to hold the comparison structures files
"""

from dataclasses import dataclass
from typing import List

from serena.utilities.ensemble_structures import Sara2SecondaryStructure


@dataclass
class ComparisonNucCounts():
    bound_count:float = -1
    unbound_count:float = -1
    both_count:float = -1
    dot_count:float = -1
    num_nucs:int = -1

@dataclass
class ComparisonResult():
    unbound_struct:Sara2SecondaryStructure 
    bound_struct: Sara2SecondaryStructure
    reference_struct: Sara2SecondaryStructure
    comp_struct: Sara2SecondaryStructure
    comp_counts: ComparisonNucCounts


class ComparisonStructures():

    def __init__(self) -> None:
        pass

    def compair_structures(self, unbound_struct:Sara2SecondaryStructure, bound_struct:Sara2SecondaryStructure, reference_struct:Sara2SecondaryStructure, nuc_count:int):
        """
        Compaire the weighted structure against the folded and not-folded mfe's.
        If a element is present in the folded mfe then it gets a '-'
        if element is in unbound only then it gets a '|'.
        The idea is that if you have a straight line in the list then it is very close to the
        folded mfe and if it is not straight then it is more like the unbound mfe.
        """
        unbound:str = '|'
        num_unbound:int = 0
        bound:str = '-'
        num_bound:int = 0
        both:str = '+'
        num_both:int = 0
        dot:str = '.'
        num_dot:int = 0
        temp_compared_struct:str = ''            

        for nuc_index in range(nuc_count):
            reference_nuc:str = reference_struct.structure[nuc_index]
            unbound_nuc:str = unbound_struct.structure[nuc_index]
            bound_nuc: str = bound_struct.structure[nuc_index]

            comp_nuc_symbol:str = ''

            if reference_nuc == bound_nuc and reference_nuc != unbound_nuc:
                comp_nuc_symbol = bound
                num_bound += 1
            elif reference_nuc != bound_nuc and reference_nuc == unbound_nuc:
                comp_nuc_symbol = unbound
                num_unbound += 1
            elif reference_nuc == bound_nuc and reference_nuc == unbound_nuc:
                comp_nuc_symbol = both
                num_both += 1
            else:
                comp_nuc_symbol = dot
                num_dot += 1
            
            temp_compared_struct = temp_compared_struct + comp_nuc_symbol
        
        comp_struct:Sara2SecondaryStructure = Sara2SecondaryStructure(sequence=unbound_struct.sequence,
                                                                        structure=temp_compared_struct)
        
        comp_nuc_counts: ComparisonNucCounts = ComparisonNucCounts(bound_count=num_bound,
                                                                    unbound_count=num_unbound,
                                                                    both_count=num_both,
                                                                    dot_count=num_dot,
                                                                    num_nucs=nuc_count)
                
        compared_result: ComparisonResult = ComparisonResult(comp_struct=comp_struct,
                                                                unbound_struct=unbound_struct,
                                                                bound_struct=bound_struct,
                                                                reference_struct=reference_struct,
                                                                comp_counts=comp_nuc_counts)
        
        return compared_result
    
class ComparisonNucAnalysis():

    def __init__(self) -> None:
        pass

    def evaluate_switchability_ensemble(self, nucs_list:NucleotideLists, settings: SwitchabilitySettings):
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
        

        for group_index in range(nucs_list.total_groups):
            bound: int = nucs_list.bound_nucs[group_index]
            unbound: int= nucs_list.unbound_nucs[group_index]
            both_nuc:int = nucs_list.both_nucs[group_index]
            dot_nuc:int = nucs_list.neithernucs[group_index]
            
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

            unbound_to_total_ratio = unbound/nucs_list.total_nucs
            bound_to_total_ratio = bound/nucs_list.total_nucs
            both_nuc_total= both_nuc/nucs_list.total_nucs
            dot_nuc_total= dot_nuc/nucs_list.total_nucs

            bound_total_list.append(bound_to_total_ratio)
            unbound_total_list.append(unbound_to_total_ratio)  

            bound_stats: str = f'BURatio:{round(bound_ratio,2)},both_Raise:{round(last_both_ratio,2)} BRaise:{round(last_bound_ratio,2)}, UDrop:{round(last_unbound_ratio,2)},BothTotal:{round(both_nuc_total,2)}, BoundTotal:{round(bound_to_total_ratio,2)}, UTotal:{round(unbound_to_total_ratio,2)}, bound_both:{round(bound_to_both_ratio,2)} B:{bound}, U:{unbound}. both:{both_nuc}'

            limit:float = settings.limit

            last_unbound_ratio = round(last_unbound_ratio,2)
            last_bound_ratio = round(last_bound_ratio,2)
            unbound_to_total_ratio = round(unbound_to_total_ratio,2)
            bound_ratio = round(bound_ratio,2)

            ev_weight_asserted:bool = nucs_list.lmv_assertions[group_index].comp_pronounced
            ev_weigth_under_limit:bool = False
            ev_weight_limit:int = 25
            if nucs_list.lmv_values[group_index].local_minima_variation_weighted_struct.ev_normalized < ev_weight_limit:
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

