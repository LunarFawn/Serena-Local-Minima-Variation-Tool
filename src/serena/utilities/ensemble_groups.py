from typing import List, Dict, NamedTuple
import struct
import pandas as pd
import sys
import openpyxl
from copy import deepcopy
from dataclasses import dataclass
from datetime import datetime, timedelta
import threading
import time
from collections import namedtuple

from serena.utilities.ensemble_structures import Sara2SecondaryStructure, Sara2StructureList, KcalRanges

class SingleEnsembleGroup():
    
    def __init__(self) -> None:
        self._group: Sara2StructureList = Sara2StructureList()
        self._multi_state_mfe_struct: List[str] = []
        """
        0 is mfe for unbound and 1 is mfe for bound
        """
        self._multi_state_mfe_kcal: List[float] = [] 
        self._kcal_span: float = 0
        self._kcal_start: float = 0
        self._kcal_end: float = 0

    @property
    def group(self):
        return self._group 

    @group.setter
    def group(self, the_group:Sara2StructureList):
        self._group = the_group
    
    @property
    def multi_state_mfe_struct(self):
        return self._multi_state_mfe_struct 

    @multi_state_mfe_struct.setter
    def multi_state_mfe_struct(self, structs:List[str]):
        self._multi_state_mfe_struct = structs
    
    def append_multi_state_mfe_data(self, structure: str, kcal: float):
        self._multi_state_mfe_struct.append(structure)
        self._multi_state_mfe_kcal.append(kcal)

    @property
    def multi_state_mfe_kcal(self):
        return self._multi_state_mfe_kcal 

    @multi_state_mfe_kcal.setter
    def multi_state_mfe_kcal(self, kcals:List[float]):
        self._multi_state_mfe_kcal = kcals

    @property
    def kcal_span(self):
        return self._kcal_span 

    @kcal_span.setter
    def kcal_span(self, kcal:float):
        self._kcal_span = kcal
    
    @property
    def kcal_start(self):
        return self._kcal_start 

    @kcal_start.setter
    def kcal_start(self, kcal:float):
        self._kcal_start = kcal

    @property
    def kcal_end(self):
        return self._kcal_end 

    @kcal_end.setter
    def kcal_end(self, kcal:float):
        self._kcal_end = kcal
    
    def update_kcals(self, start:float, stop:float, span:float):
        self._kcal_start = start
        self._kcal_end = stop
        self._kcal_span = span

class MultipleEnsembleGroups():

    def __init__(self, non_switch_kcal:float =0, non_switch_struct:Sara2SecondaryStructure= Sara2SecondaryStructure(), switched_kcal:float=0, switched_struct:Sara2SecondaryStructure=Sara2SecondaryStructure()) -> None:
        self._groups: List[SingleEnsembleGroup] = []  
        self._raw_groups: List[Sara2StructureList] = []
        self._non_switch_state_mfe_kcal: float = non_switch_kcal
        self._non_switch_state_structure: Sara2SecondaryStructure = non_switch_struct
        self._switched_state_mfe_kcal: float = switched_kcal
        self._switched_state_structure: Sara2SecondaryStructure = switched_struct
        self._groups_dict: Dict[int, Sara2StructureList] = {}
        self._group_values: List[float] = []
        self._num_groups: int = 0
        self._group_kcal_ranges: List[KcalRanges] =  []
    
    @property
    def num_groups(self):
        return self._num_groups
    
    @num_groups.setter
    def num_groups(self, num: int):
        self._num_groups = num

    #def add_group(self, group:SingleEnsembleGroup, group_index:int, value_of_group:float, start_kcal:float = 0, end_kcal:float=0):
    #    if self._switched_state_mfe_kcal >= group.kcal_start and self._switched_state_mfe_kcal < group.kcal_end:
    #        group.has_bound_mfe_kcal = True
    #    self._groups.append(group)
    #    self._raw_groups.append(group.group)
    #    self._groups_dict[group_index]= group.group
    #    self._group_values.append(value_of_group)
    #    kcal_range: KcalRanges = KcalRanges(start=start_kcal, stop=end_kcal)
    #    self._group_kcal_ranges.append(kcal_range)
    
    #def append_group(self, group:SingleEnsembleGroup, group_value: float, start_kcal:float = 0, end_kcal:float=0):
    #    self._num_groups = self._num_groups + 1
    #    self._groups.append(group)
    #    self._raw_groups.append(group.group)
    #    self._groups_dict[self._num_groups-1]= group.group
    #   self._group_values.append(group_value)
    #    kcal_range: KcalRanges = KcalRanges(start=start_kcal, stop=end_kcal)
    #    self._group_kcal_ranges.append(kcal_range)

    @property
    def groups(self):
        return self._groups
    
    @groups.setter
    def groups(self, groupss:List[SingleEnsembleGroup]):
        self._groups = groupss
    
    @property
    def raw_groups(self):
        return self._raw_groups
    
    @raw_groups.setter
    def raw_groups(self, groupss:List[Sara2StructureList]):
        self._raw_groups = groupss
    
    @property
    def non_switch_state_mfe_kcal(self):
        return self._non_switch_state_mfe_kcal
    
    @property
    def non_switch_state_structure(self)->Sara2SecondaryStructure:
        return self._non_switch_state_structure
    
    @property
    def switched_state_mfe_kcal(self):
        return self._switched_state_mfe_kcal
    
    @property
    def switched_state_structure(self)->Sara2SecondaryStructure:
        return self._switched_state_structure
    
    @property
    def groups_dict(self):
        return self._groups_dict
    
    @groups_dict.setter
    def groups_dict(self, dict: Dict[int, Sara2StructureList]):
        self._groups_dict = dict
    
    @property
    def group_values(self):
        return self._group_values
    
    @group_values.setter
    def group_values(self, values :List[float]):
        self._group_values = values
    
    @property
    def group_kcal_ranges(self):
        return self._group_kcal_ranges
    
    @group_kcal_ranges.setter
    def group_kcal_ranges(self, values :List[KcalRanges]):
        self._group_kcal_ranges = values

    @property
    def total_structures(self):
        total:int = 0
        for group in self.raw_groups:
            total += group.num_structures
        return total
    


