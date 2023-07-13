

from dataclasses import dataclass

@dataclass
class SwitchAnalysisDataToProcess():
    sequence:str = ''
    span: int = -1
    units:float = -1
    fmn_struct:str = ''
    fmn_struct_free_energy:float = -1
    temperature:int = -1