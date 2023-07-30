"""
File to handles teh calsses for dealing with scores
"""

from dataclasses import dataclass
from typing import List

from serena.utilities.judge_pool import JudgesResults
from serena.utilities.investigator import InvestigatorResults

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

    def score_groups(self, judge_results:JudgesResults, investigator: InvestigatorResults):
        #inititalization data
        found_functional_switch: List[int] = judge_results.switchable_groups_list 
        found_powerful_switch: List[int] = judge_results.powerfull_groups_list
        found_on_off_switch: List[int] = judge_results.on_off_groups_list


        bound_range_index_plus_one:List[int] = judge_results.switchable_groups_list
        is_powerful_switch:bool = judge_results.is_powerful_switch
        is_functional_switch:bool = judge_results.is_good_switch
        is_off_on_switch:bool = judge_results.is_on_off_switch

        #SetupScores
        total_score:float = 0
        functional_switch_score:float = 0
        powerful_switch_score:float = 0
        on_off_switch_score:float = 0
        bonuses:float = 0
        penalties:float = 0

        #main scores
        if is_powerful_switch is True:
            multiplier:int = 1
            message:str = 'Potential High Fold Change'
            #result_messages = self.log_message(message, result_messages) 
            powerful_switch_score = powerful_switch_score + (len(found_powerful_switch) * multiplier)
        
        if is_functional_switch is True: 
            multiplier:int = 1
            message:str = "Potential  Functional Switch"
            #result_messages = self.log_message(message, result_messages)
            functional_switch_score = functional_switch_score + (len(found_functional_switch) * multiplier)
        
        if is_off_on_switch is True:
            multiplier:int = 1
            message:str = "Potential  off/on leaning design via LMV"
            #result_messages = self.log_message(message, result_messages)
            on_off_switch_score= on_off_switch_score + (len(found_on_off_switch) * multiplier)

        total_score = powerful_switch_score + functional_switch_score + on_off_switch_score

        #now bonuses
        for value in found_functional_switch:
            if value >= 0 and value <= 1 and value != -1:
                message:str = "Confirmned good. Add bonus point for point for functional being in first two groups"
                #result_messages = self.log_message(message, result_messages)
                functional_switch_score += 1
                bonuses += 1

            if value in found_on_off_switch:
                message:str = "Add bonus for functional being in range of on/off prediction"
                #result_messages = self.log_message(message, result_messages)
                functional_switch_score += 1
                bonuses += 1

        for value in found_powerful_switch:
            if value >= 0 and value <= 1 and value != -1:
                message:str = "Confirmned good. Add bonus point for high performing being in first two groups"
                #result_messages = self.log_message(message, result_messages)
                powerful_switch_score += 1
                bonuses += 1

            if value in found_on_off_switch:
                message:str = "Add bonus for high performing being in range of on/off prediction"
                #result_messages = self.log_message(message, result_messages)
                powerful_switch_score += 1
                bonuses += 1
      
        #not sure if I want to use... i was ify about before and it seams not fully baked in implementatiuon. 
        # need to make a ticket for this funciton
        #if is_good_switch is True and bound_to_both_ratio >= 0.08:
        #    message:str = "Low number of both and mfe nucs in relation to bound. Add bonus point"
        #    result_messages = self.log_message(message, result_messages)
        #    score= score + 1

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
