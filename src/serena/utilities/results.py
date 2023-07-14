
from dataclasses import dataclass, fields, Field
from typing import List
import pandas as pd 
from pandas import DataFrame

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

class ResultsToDataframe():
    
    def __init__(self) -> None:
        pass

    def make_switch_analysis_dataframe(self, result: SwitchAnalysisResults):
        header:List[str] = []
        values:List[str] = []
        for field in fields(result):
            field_name = field.name
            header.append(field_name)
            filed_value = getattr(result, field_name)
            values.append(str(filed_value))
        df_result:DataFrame = pd.DataFrame(values, columns=header)
        return df_result

