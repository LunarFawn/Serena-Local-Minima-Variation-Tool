
from dataclasses import dataclass


@dataclass
class SwitchAnalysisResults():
    bound_over_unbound_ratio:float = -1
    both_last_raise_ratio:float = -1 
    bound_last_raise_ratio:float = -1 
    unbound_last_drop_ratio:float = -1
    both_over_total_ratio:float = -1
    bound_over_total_ratio: float = -1
    unbound_over_total_ratio: float = -1
    bound_over_both_minus_unbound_ratio:float = -1
    bound_total:int = -1
    unbound_total:int = -1
    both_total:int =-1
    ev_comparison: float = -1
    ev_relative: float = -1
    ev_unbound_mfe:float = -1

