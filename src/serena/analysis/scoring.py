"""
File to handles teh calsses for dealing with scores
"""

from dataclasses import dataclass
from typing import List

from serena.analysis.judge_pool import JudgesResults
from serena.analysis.investigator import InvestigatorResults

@dataclass
class BasicScoreResults():
    """
    Basic scores for switchyness
    """
    total_score:float = 0
    functional_switch_score:float = 0
    powerful_switch_score:float = 0
    on_off_switch_score:float = 0
    bonuses:float = 0
    penalties:float = 0

@dataclass
class AdvancedScoreResults():
    """
    Bonuses that amplify ability to decern switchyness
    """
    lmv_bonus:float =0
    lmv_penalty:float = 0    
    comp_bonus:float = 0
    comp_penalty:float = 0
    excess_struct_penalty:float=0
    total_score:float = 0


class SerenaScoring():
    """
    Scoring the results from the judges decisions
    """
    def __init__(self) -> None:
        pass

    def basic_score_groups(self, judge_results:JudgesResults, investigator: InvestigatorResults)->BasicScoreResults:
        """
        Perform basic scoring functions that determine switchyness of rna sequence
        """
        #inititalization data
        found_functional_switch: List[int] = judge_results.comp_switch_judge.switchable_groups_list 
        found_powerful_switch: List[int] = judge_results.comp_switch_judge.powerfull_groups_list
        found_on_off_switch: List[int] = judge_results.lmv_switch_judge.on_off_groups_list


        bound_range_index_plus_one:List[int] = judge_results.comp_switch_judge.switchable_groups_list
        is_powerful_switch:bool = judge_results.comp_switch_judge.is_powerful_switch
        is_functional_switch:bool = judge_results.comp_switch_judge.is_good_switch
        is_off_on_switch:bool = judge_results.lmv_switch_judge.is_on_off_switch

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
        
        total_score = powerful_switch_score + functional_switch_score + on_off_switch_score - penalties

        basic_score_results:BasicScoreResults = BasicScoreResults(total_score=total_score,
                                                                  functional_switch_score=functional_switch_score,
                                                                  powerful_switch_score=powerful_switch_score,
                                                                  on_off_switch_score=on_off_switch_score,
                                                                  bonuses=bonuses,
                                                                  penalties=penalties)

        return basic_score_results
      
    def excessive_structures_penalties(self, num_structures: int, excess_divisor:float,excess_limit:float):
        """
        Algorithm for determining the penalty for excessive number of secondary structures in the 
        whole ensemble
        """
        #excess_divisor:float = 2000#2500
        penalty:float = 0
        if num_structures > excess_limit:
            factor:float = ((float(num_structures) - excess_limit) / excess_divisor ) * .5
            message:str = f'Exsessive structs. Found:{num_structures} penalizing {factor} points '
            #result_messages = self.log_message(message, result_messages)
            sixty_range_num:float = 50000#15000
            #penalize for too many structs
            penalty += factor
            if num_structures > sixty_range_num:
                message:str = f'Significant excess structures found: found {num_structures - sixty_range_num} structures over limit of {sixty_range_num}'
                #result_messages = self.log_message(message, result_messages)
                message:str = f'Eterna_score should be ~60 for temp group and could be good design currently has high penalty for excess structures and now yet one more penalty'
                #result_messages = self.log_message(message, result_messages)
                penalty += .5
        
        return penalty
    
    def advanced_score_groups(self, judge_results:JudgesResults, investigator: InvestigatorResults):
        """
        Bonuses and penalties that affect the fine tunning of swithyness determinations
        """
        lmv_bonus:float =0
        lmv_penalty:float = 0    
        comp_bonus:float = 0
        comp_penalty:float = 0
        total_score:float = 0

        comp_less_ratio: float = investigator.lmv_assertions.comp_compare_to_mfe.count('<') / investigator.num_groups
        com_great_ratio: float = investigator.lmv_assertions.comp_compare_to_mfe.count('>')  / investigator.num_groups
        message:str = f'ev comp great:{com_great_ratio}, ev comp less:{comp_less_ratio}'
        #result_messages = self.log_message(message, result_messages)
        if com_great_ratio < comp_less_ratio and comp_less_ratio >= .7:
            message:str = "EV for comparison struct is LESS MORE OFTEN than unbound mfe so add bonus"
            #result_messages = self.log_message(message, result_messages)
            lmv_bonus += 1
        elif com_great_ratio > comp_less_ratio and com_great_ratio >= .5:
            message:str = "EV for comparison struct is GREATER MORE OFTEN than unbound mfe so penatly"
            #result_messages = self.log_message(message, result_messages)
            lmv_penalty += .5
            if com_great_ratio >= .8:
                message:str = "EV for comp is GREATER EXTRA MORE OFTEN then mfe so minus penalty point"
                #result_messages = self.log_message(message, result_messages)
                lmv_penalty += .5
        
        if investigator.comparison_eval_results.nuc_penatly_count > 0:
            if investigator.comparison_eval_results.BUratio_list[0] >= .75:
                new_penalty: float = investigator.comparison_eval_results.nuc_penatly_count * .5
                message:str = f'Bound unbound ratio higher than 75% so it will most likely just fold into what should have been a switch so minus {new_penalty} points'
                #result_messages = self.log_message(message, result_messages)
                comp_penalty += new_penalty
            #elif BUratio_list[0] > .60 and BUratio_list[1] < .3:
            #    new_penalty: float = nuc_penatly_count * 1
            #    message:str = f'Bound unbound ratio higher than 50% and then the 2nd energy group less than 20% so it will likely be blocked from switching so minus {new_penalty} points'
            #    result_messages = self.log_message(message, result_messages)
            #    score = score - new_penalty
            else:
                new_penalty: float = investigator.comparison_eval_results.nuc_penatly_count * .5                   
                message:str = f'Bound nucs found in first energy group. Design is primed to switch so add bonus of {new_penalty} points'
                #result_messages = self.log_message(message, result_messages)
                comp_bonus += new_penalty

        #not sure if I want to use... i was ify about before and it seams not fully baked in implementatiuon. 
        # need to make a ticket for this funciton
        #if is_good_switch is True and bound_to_both_ratio >= 0.08:
        #    message:str = "Low number of both and mfe nucs in relation to bound. Add bonus point"
        #    result_messages = self.log_message(message, result_messages)
        #    score= score + 1
        excess_limit:float = 20000#this is based on new data 7500
        excess_divisor:float = 2000#2500
        excess_struct_penalty:float = self.excessive_structures_penalties(num_structures=investigator.total_structures_ensemble,
                                                                    excess_limit=excess_limit,
                                                                    excess_divisor=excess_divisor)
        total_score = lmv_bonus - lmv_penalty + comp_bonus - comp_penalty - excess_struct_penalty
        advanced_score_response: AdvancedScoreResults = AdvancedScoreResults(lmv_bonus=lmv_bonus,
                                                                             lmv_penalty=lmv_penalty,
                                                                             comp_bonus=comp_bonus,
                                                                             comp_penalty=comp_penalty,
                                                                             excess_struct_penalty=excess_struct_penalty,
                                                                             total_score=total_score)
        return advanced_score_response

    