"""
Sara2 api for accessing and manipulating secondary structures 
in dot parenthisis form
copyright 2023 GrizzlyEngineer
"""
from typing import List, Dict
import struct
import pandas as pd
import sys
import openpyxl
from copy import deepcopy
from dataclasses import dataclass
from datetime import datetime, timedelta
import threading
import time




class Sara2SecondaryStructure(object):

    def __init__(self, sequence:str = '', structure: str = '', freeEnergy: float = 0, stackEnergy: float = 0) -> None:
        self._sequence: str = sequence
        self._structure: str = structure
        self._freeEnergy: float = freeEnergy
        self._stackEnergy: float = stackEnergy
        #self._nuc_count: int = len(sequence)

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, primary_struc: str):
        self._sequence = primary_struc

    @property
    def structure(self):
        return self._structure

    @structure.setter
    def structure(self, dot_parens: str):
        self._structure = dot_parens
    
    @property
    def freeEnergy(self):
        return self._freeEnergy

    @freeEnergy.setter
    def freeEnergy(self, energy: float):
        self._freeEnergy = energy
    
    @property
    def stackEnergy(self):
        return self._stackEnergy

    @stackEnergy.setter
    def stackEnergy(self, energy: float):
        self._stackEnergy = energy
    
    @property
    def nuc_count(self):
        return len(self._sequence)

class ComparisonStructures():

    def __init__(self) -> None:
        self._structures: Sara2StructureList = Sara2StructureList()
        self._names: List[str] = []
        self._structures_dict: Dict[str, Sara2SecondaryStructure]

    def add_structure(self, structure:Sara2SecondaryStructure, name: str):
        self._structures.add_structure(structure)
        self._names.append(name)
        self._structures_dict[name] = structure
    
    @property
    def structures(self):
        return self._structures
    
    @property
    def names(self):
        return self._names
    
    @property
    def structures_dict(self):
        return self._structures_dict
    
    def get_structure_by_name(self, name:str):
        structure: Sara2SecondaryStructure = Sara2SecondaryStructure()
        if name in self._structures_dict:
            structure = self._structures_dict[name]
        else:
            raise Exception("name does not exist")       
        return structure
    
    def get_structure_by_index(self, index:int):
        structure: Sara2SecondaryStructure = Sara2SecondaryStructure()
        
        if index <= len(self._structures)-1:
            structure = self._structures[index]
        else:
            raise Exception("wrong index")
                
        return structure
    

class Sara2StructureList(object):
    
    def __init__(self) -> None:
        self._sara_structures_list: List[Sara2SecondaryStructure] = []
        self._structures: List[str] = []
        self._freeEnergy_list: list[float] = []
        self._stackEnergy_list: list[float] = []
        self._min_freeEnergy: float = 0
        self._max_freeEnergy: float = 0
        self._min_stackEnergy: float = 0
        self._max_stackEnergy: float = 0
        self._num_structures: int = 0
        self._nuc_count: int = 0
        self._mfe_structure: str = ''
        self._mfe_freeEnergy: float = 0
        self._mfe_stackEnergy: float = 0
        self._freeEnergy_span:float = 0
        self._stackEnergy_span:float = 0
    

    def process_energy(self):
            #now populate min and max
        #do free energy
        if len(self._freeEnergy_list) == 0:
            self._min_freeEnergy = 0
            self._max_freeEnergy = 0
        else:
            self._min_freeEnergy = min(self._freeEnergy_list)
            self._max_freeEnergy = max(self._freeEnergy_list)

        self._freeEnergy_span = self._max_freeEnergy - self._min_freeEnergy
        #do stack energy

        if len(self._stackEnergy_list) == 0:
            self._min_stackEnergy = 0
            self._max_stackEnergy = 0
        else:
            self._min_stackEnergy = min(self._stackEnergy_list)
            self._max_stackEnergy = max(self._stackEnergy_list)
        self._stackEnergy_span = self._max_stackEnergy - self._min_stackEnergy

        #now count
        if len(self._structures) == 0:
            self._num_structures = 0
        else:
            self._num_structures = len(self._structures)

    def add_structure(self, structure: Sara2SecondaryStructure):
        self._sara_structures_list.append(structure)
        self._structures.append(structure.structure)
        self._freeEnergy_list.append(structure.freeEnergy)
        self._stackEnergy_list.append(structure.stackEnergy)
        #self.process_energy()


    def remove_structure(self, index:int):
        del self._structures[index]
        del self._freeEnergy_list[index]
        del self._stackEnergy_list[index]
        #self.process_energy()            


    @property
    def mfe_structure(self):
        return self.sara_stuctures[0].sequence

    @property
    def mfe_freeEnergy(self):
        return self.sara_stuctures[0].freeEnergy
    
    @property
    def mfe_stackEnergy(self):
        return self.sara_stuctures[0].stackEnergy
    
    @property
    def nuc_count(self):
        return self.sara_stuctures[0].nuc_count 

    @property
    def sara_stuctures(self):
        return self._sara_structures_list

    @sara_stuctures.setter
    def sara_stuctures(self, structs_list: List[Sara2SecondaryStructure]):
        #reset list
        self._sara_structures_list=[]
        #fill it in now
        for struc in structs_list:
            self.add_structure(struc)

    
    @property
    def max_free_energy(self):
        self.process_energy()
        return self._max_freeEnergy
    
    @property
    def min_free_energy(self):
        self.process_energy()
        return self._min_freeEnergy
    
    @property
    def max_stack_energy(self):
        self.process_energy()
        return self._max_stackEnergy
    
    @property
    def min_stack_energy(self):
        self.process_energy()
        return self._min_stackEnergy
    
    @property
    def num_structures(self):
        self.process_energy()
        return self._num_structures
    
    @property
    def freeEnergy_span(self):
        self.process_energy()
        return self._freeEnergy_span

    @property
    def stackEnergy_span(self):
        self.process_energy()
        return self._stackEnergy_span 
    
    @property
    def weighted_structure(self):
        return self._weighted_structure
    
    @weighted_structure.setter
    def weighted_structure(self, structure: str):
        self._weighted_structure = structure


@dataclass
class EV:
    ev_normalized: float = -1
    ev_ThresholdNorm: float = -1
    ev_structure: float = -1

@dataclass
class EVResult():
    groups_list : List[Sara2StructureList]
    groups_dict: Dict[int, Sara2StructureList]
    group_values: List[float]
    group_ev_list: List[EV]
    group_ev_dict: Dict[int,EV]


class LMV_Token():
    def __init__(self, num_groups: int) -> None:
        self._group_results: List[EV] = num_groups * [EV()]
        self._group_dict: Dict[int,EV] = {}
        self._group_values: List[str] = num_groups * ['']
        self._group_done_status: List[bool] = num_groups * [False]
    
    @property
    def group_dict(self):
        return self._group_dict
        
    def set_group_dict(self, index:int, value:EV):
        self._group_dict[index]=value

    @property
    def group_results(self):
        return self._group_results
        
    def set_group_result(self, index:int, value:EV):
        self._group_results[index]=value
    
    @property
    def group_values(self):
        return self._group_values
        
    def set_group_values(self, index:int, value:str):
        self._group_values[index]=value

    @property
    def group_done_status(self):
        return self._group_done_status
        
    def set_group_done_status(self, index:int, state:bool):
        self._group_done_status[index]=state
    
    @property
    def is_done(self):
        is_completed:bool = False
        if self._group_done_status.count(False) == 0:
            #its done
            is_completed = True
        return is_completed


class LMV_Shuttle():

    def __init__(self, structs_list: Sara2StructureList, mfe:Sara2SecondaryStructure, group_index:int, token:LMV_Token) -> None:
        self._kcal_group_structures_list: Sara2StructureList = structs_list
        self._sara_mfestructure:Sara2SecondaryStructure = mfe
        self._group_index:int = group_index
        self._token:LMV_Token = token
    
    @property
    def kcal_group_structures_list(self):
        return self._kcal_group_structures_list

    @kcal_group_structures_list.setter
    def kcal_group_structures_list(self, new_list: Sara2StructureList):
        self._kcal_group_structures_list = new_list
    
    @property
    def sara_mfestructure(self):
        return self._sara_mfestructure

    @sara_mfestructure.setter
    def sara_mfestructure(self, new_strucr: Sara2SecondaryStructure):
        self._sara_mfestructure = new_strucr
    
    @property
    def group_index(self):
        return self._group_index

    @group_index.setter
    def group_index(self, new_index: int):
        self._group_index = new_index
    
    @property
    def token(self):
        return self._token

    @token.setter
    def token(self, new_token: LMV_Token):
        self._token = new_token


#@dataclass
class SingleEnsembleGroup():
    
    def __init__(self) -> None:
        self._group: Sara2StructureList = Sara2StructureList()
        self._multi_state_mfe_struct: List[str] = []
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
        return self.kcal_span 

    @kcal_span.setter
    def kcal_span(self, kcal:float):
        self._kcal_span = kcal
    
    @property
    def kcal_start(self):
        return self.kcal_start 

    @kcal_start.setter
    def kcal_start(self, kcal:float):
        self._kcal_start = kcal

    @property
    def kcal_end(self):
        return self.kcal_end 

    @kcal_end.setter
    def kcal_end(self, kcal:float):
        self._kcal_end = kcal
    
    def update_kcals(self, start:float, stop:float, span:float):
        self._kcal_start = start
        self._kcal_end = stop
        self._kcal_span = span

class MultipleEnsembleGroups():

    def __init__(self, non_switch_kcal:float =0, non_switch_struct:str = '', switched_kcal:float=0, switched_struct:str='') -> None:
        self._groups: List[SingleEnsembleGroup] = []  
        self._raw_groups: List[Sara2StructureList] = []
        self._non_switch_state_mfe_kcal: float = non_switch_kcal
        self._non_switch_state_structure: str = non_switch_struct
        self._switched_state_mfe_kcal: float = switched_kcal
        self._switched_state_structure: str = switched_struct
        self._groups_dict: Dict[int, Sara2StructureList] = {}
        self._group_values: List[float] = []
        self._num_groups: int = 0
    
    @property
    def num_groups(self):
        return self._num_groups
    
    @num_groups.setter
    def num_groups(self, num: int):
        self._num_groups = num

    def add_group(self, group:SingleEnsembleGroup, group_index:int, value_of_group:float):
        if self._switched_state_mfe_kcal >= group.kcal_start and self._switched_state_mfe_kcal < group.kcal_end:
            group.has_bound_mfe_kcal = True
        self._groups.append(group)
        self._raw_groups.append(group.group)
        self._groups_dict[group_index]= group.group
        self._group_values.append(value_of_group)
    
    def append_group(self, group:SingleEnsembleGroup, group_value: float):
        self._num_groups = self._num_groups + 1
        self._groups.append(group)
        self._raw_groups.append(group.group)
        self._groups_dict[self._num_groups-1]= group.group
        self._group_values.append(group_value)

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
    def non_switch_state_structure(self):
        return self._non_switch_state_structure
    
    @property
    def switched_state_mfe_kcal(self):
        return self._switched_state_mfe_kcal
    
    @property
    def switched_state_structure(self):
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

@dataclass
class WeightedStructureData():
    raw_group: Sara2StructureList = Sara2StructureList()
    weighted_dot_paren_structure: str = ''
    weighted_compared_line:str = ''
    unbound_mfe_dot_paren_struct: str = ''
    unbound_mfe_kcal:float = 0
    bound_mfe_dot_paren_struct: str = ''
    bound_mfe_kcal:float = 0
    BURatio: float = -1
    BRaise: float = -1
    UDrop: float = -1
    UTotal: float = -1
    bound_num:float = -1
    unbound_num: float = -1
    switch_score:float = -1
    kcal_start:float = -1
    kcal_stop:float = -1
