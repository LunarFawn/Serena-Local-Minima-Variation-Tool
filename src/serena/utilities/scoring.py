"""
File to handles teh calsses for dealing with scores
"""

from dataclasses import dataclass
from typing import List

from serena.utilities.judge_pool import JudgesResults

class SpecialPenalties():

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

class SerenaScoring():
    def __init__(self) -> None:
        pass

    def score_groups(self, judge_results:JudgesResults):            
        bound_range_index_plus_one:List[int] = judge_results.switchable_groups_list
        is_powerful_switch:bool = judge_results.is_powerful_switch
        is_good_switch:bool = judge_results.is_good_switch

        bound_range_min_minus_1: int = 0
        bound_range_max_plus: int = 0
        bound_good_range_modifier:int = 0 
        #bound_range_index_plus_one is the list for groups with pair
        if len(bound_range_index_plus_one) > 0: 
            bound_range_min_minus_1: int = min(bound_range_index_plus_one) - bound_good_range_modifier
        if len(bound_range_index_plus_one) > 0: 
            bound_range_max_plus: int = max(bound_range_index_plus_one) + bound_good_range_modifier            

        if is_powerful_switch is True:
            message:str = 'Potential High Fold Change'
            result_messages = self.log_message(message, result_messages) 
            score = score + 1
        
        if is_good_switch is True: 
            message:str = "Potential  Functional Switch"
            result_messages = self.log_message(message, result_messages)
            score = score + (len(found_bound_ratio_list)*1)
        
        if is_off_on_switch is True:
            message:str = "Potential  off/on leaning design via LMV"
            result_messages = self.log_message(message, result_messages)
            score= score + 1
        
        #this is lmv stuff
        if found_bound_index >= bound_range_min_minus_1 and found_bound_index <= bound_range_max_plus and found_bound_index != -1 and is_off_on_switch is True:
            message:str = "Confirmned good. Add bonus point for on/off via LMV being in range for folding"
            result_messages = self.log_message(message, result_messages)
            score= score + 1
        elif found_bound_index <= 2 and found_bound_index != -1 and is_in_bound_range is True:
            message:str = "Confirmned good. Add bonus point for on/off via LMV being in first three groups"
            result_messages = self.log_message(message, result_messages)
            score= score + 1
        for value in found_bound_ratio_list:
            if value >= bound_range_min_minus_1 and value <= bound_range_max_plus and found_bound_ratio_index != -1:
                message:str = "Confirmned good. Add bonus point for functional being in range for folding"
                result_messages = self.log_message(message, result_messages)
                score= score + 1
            elif value >= 0 and value <= 1 and value != -1:
                message:str = "Confirmned good. Add bonus point for point for functional being in first two groups"
                result_messages = self.log_message(message, result_messages)
                score= score + 1

        if found_bound_ratio_high_index >= bound_range_min_minus_1 and found_bound_ratio_high_index <= bound_range_max_plus and found_bound_ratio_high_index != -1 :
            message:str = "Confirmned good. Add bonus point for high performing being in range for folding"
            result_messages = self.log_message(message, result_messages)
            score= score + 1
        elif found_bound_ratio_high_index >= 0 and found_bound_ratio_high_index <= 1 and found_bound_ratio_high_index != -1:
            message:str = "Confirmned good. Add bonus point for high performing being in first two groups"
            result_messages = self.log_message(message, result_messages)
            score= score + 1

        if found_bound_ratio_high_index in found_bound_list:
            message:str = "Add bonus for high performing being in range of on/off prediction"
            result_messages = self.log_message(message, result_messages)
            score= score + 1
        
        if found_bound_ratio_index in found_bound_list:
            message:str = "Add bonus for functional being in range of on/off prediction"
            result_messages = self.log_message(message, result_messages)
            score= score + 1

        excess_limit:float = 20000#this is based on new data 7500
        if span_structures.num_structures > excess_limit:#15000:
            excess_divisor:float = 2000#2500
            factor:float = ((float(span_structures.num_structures) - excess_limit) / excess_divisor ) * .5
            message:str = f'Exsessive structs. Found:{span_structures.num_structures} penalizing {factor} points '
            result_messages = self.log_message(message, result_messages)
            sixty_range_num:float = 50000#15000
            #penalize for too many structs
            score = score - factor
            if span_structures.num_structures > sixty_range_num:
                message:str = f'Significant excess structures found: found {span_structures.num_structures - sixty_range_num} structures over limit of {sixty_range_num}'
                result_messages = self.log_message(message, result_messages)
                message:str = f'Eterna_score should be ~60 for temp group and could be good design currently has high penalty for excess structures and now yet one more penalty'
                result_messages = self.log_message(message, result_messages)
                score = score - .5
        
        if is_good_switch is True and bound_to_both_ratio >= 0.08:
            message:str = "Low number of both and mfe nucs in relation to bound. Add bonus point"
            result_messages = self.log_message(message, result_messages)
            score= score + 1

        comp_less_ratio: float = ev_comp_to_mfe.count('<') / num_groups
        com_great_ratio: float = ev_comp_to_mfe.count('>')  / num_groups
        message:str = f'ev comp great:{com_great_ratio}, ev comp less:{comp_less_ratio}'
        result_messages = self.log_message(message, result_messages)
        if com_great_ratio < comp_less_ratio and comp_less_ratio >= .7:
            message:str = "EV for comparison struct is LESS MORE OFTEN than unbound mfe so add bonus"
            result_messages = self.log_message(message, result_messages)
            score= score + 1
        elif com_great_ratio > comp_less_ratio and com_great_ratio >= .5:
            message:str = "EV for comparison struct is GREATER MORE OFTEN than unbound mfe so penatly"
            result_messages = self.log_message(message, result_messages)
            score= score - .5
            if com_great_ratio >= .8:
                message:str = "EV for comp is GREATER EXTRA MORE OFTEN then mfe so minus penalty point"
                result_messages = self.log_message(message, result_messages)
                score= score - .5
        
        if nuc_penatly_count > 0:
            if BUratio_list[0] >= .75:
                new_penalty: float = nuc_penatly_count * .5
                message:str = f'Bound unbound ratio higher than 75% so it will most likely just fold into what should have been a switch so minus {new_penalty} points'
                result_messages = self.log_message(message, result_messages)
                score = score - new_penalty
            #elif BUratio_list[0] > .60 and BUratio_list[1] < .3:
            #    new_penalty: float = nuc_penatly_count * 1
            #    message:str = f'Bound unbound ratio higher than 50% and then the 2nd energy group less than 20% so it will likely be blocked from switching so minus {new_penalty} points'
            #    result_messages = self.log_message(message, result_messages)
            #    score = score - new_penalty
            else:
                new_penalty: float = nuc_penatly_count * .5                   
                message:str = f'Bound nucs found in first energy group. Design is primed to switch so add bonus of {new_penalty} points'
                result_messages = self.log_message(message, result_messages)
                score = score + new_penalty
