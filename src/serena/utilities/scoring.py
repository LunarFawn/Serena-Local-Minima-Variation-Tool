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