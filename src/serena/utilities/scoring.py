"""
File to handles teh calsses for dealing with scores
"""


class Penalties():

    def __init__(self) -> None:
        pass

    def excessive_structures(self, num_structures: int, excess_divisor:float,excess_limit:float):
        #excess_divisor:float = 2000#2500
        score:float = 0
        factor:float = ((float(num_structures) - excess_limit) / excess_divisor ) * .5
        message:str = f'Exsessive structs. Found:{num_structures} penalizing {factor} points '
        result_messages = self.log_message(message, result_messages)
        sixty_range_num:float = 50000#15000
        #penalize for too many structs
        score = score - factor
        if num_structures > sixty_range_num:
            message:str = f'Significant excess structures found: found {num_structures - sixty_range_num} structures over limit of {sixty_range_num}'
            result_messages = self.log_message(message, result_messages)
            message:str = f'Eterna_score should be ~60 for temp group and could be good design currently has high penalty for excess structures and now yet one more penalty'
            result_messages = self.log_message(message, result_messages)
            score = score - .5
        
        return score

class Bonuses():
    def __init__(self) -> None:
        pass

@dataclass
class SwitchabilitySettings():
    limit: float = 1.5 

@dataclass
class NucleotideLists():
    unbound_nucs:List[int]
    bound_nucs:List[int]
    both_nucs:List[int]
    neithernucs:List[int]
    lmv_assertions:List[LMVResult]
    lmv_values:List[WeightedLocalMinimaVariation]
    total_nucs: int
    total_groups:int

class Evaluations():

    def __init__(self) -> None:
        pass

    def evaluate_comparison_ensemble(self, nucs_list:NucleotideLists, settings: SwitchabilitySettings):
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
