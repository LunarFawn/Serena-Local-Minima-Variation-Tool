
from dataclasses import dataclass
from typing import List, Dict
from enum import Enum
import threading
from datetime import datetime

from serena.apps.original_weighted_analysis import Sara2SecondaryStructure, Sara2StructureList, EnsembleVariation, EVResult





@dataclass
class SwitchInput():
    sequence:str
    fmn_struct:str
    fmn_struct_free_energy:float
    span:int
    units:int
    run_name:str


class Temperature_Shuttle():
    pass

class Temperature_Token():
    pass

class Temperature_ThreadProcessor():
    
    def __init__(self, temp_list: List[int], switch_input:SwitchInput) -> None:
        self._temp_list: List[int] = temp_list
        num_temps:int = len(temp_list)
        self._num_temps: int =  num_temps
        self._temp_token: Temperature_Token = Temperature_Token(num_temps)
        self._SwitchAnalysis: EnsembleVariation = EnsembleVariation()
        self._switch_input: SwitchInput = switch_input

    @property
    def temp_list(self):
        return self._temp_list

    @temp_list.setter
    def temp_list(self, temp_list:List[int]):
        self._temp_list = temp_list
    
    @property
    def switch_input(self):
        return self._switch_input

    @switch_input.setter
    def switch_input(self, new_input:SwitchInput):
        self._switch_input = new_input
    
    @property
    def num_temps(self):
        return self._num_temps

    @num_temps.setter
    def num_temps(self, new_num:int):
        self._num_temps = new_num

    @property
    def temp_token(self):
        return self._temp_token

    @temp_token.setter
    def temp_token(self, new_token:Temperature_Token):
        self._temp_token = new_token
    
    @property
    def SwitchAnalysis(self):
        return self._SwitchAnalysis

    @SwitchAnalysis.setter
    def SwitchAnalysis(self, new_lmv:EnsembleVariation):
        self._SwitchAnalysis = new_lmv

    def run_SwitchAnalysis(self):
        self.start_calculations()
        self.wait_for_finish()
        #the test should be done now
        #check foor index that is -1 and if so then use prev value
        num_groups:int = len(self.temp_token.temps_results)
        #for index in range(1, num_groups):
        #    if self.temp_token.temps_results[index].ev_normalized == -1:
        #        previous_EV = self.group_token.group_results[index-1]
        #        self.group_token.group_results[index] = previous_EV
        #        self.group_token.group_dict[index] = previous_EV
        return self.temp_token

    def start_calculations(self):     
        for thread_index in range(self.num_temps):
            sequence:str = self.switch_input.sequence
            fmn_struct:str = self.switch_input.fmn_struct
            fmn_struct_free_energy:float = self.switch_input.fmn_struct_free_energy
            span:int = self.switch_input.span
            units:int = self.switch_input.units
            run_name:str = self.switch_input.run_name
            new_shuttle: Temperature_Shuttle = Temperature_Shuttle(sequence=self.switch_input.sequence,
                                                                fmn_struct=self.switch_input.fmn_struct,
                                                                fmn_struct_free_energy=self.switch_input.fmn_struct_free_energy,
                                                                span=self.switch_input.span,
                                                                units=self.switch_input.units,
                                                                group_index=thread_index,
                                                                token=self.temp_token) 
            mew_thread = threading.Thread(target=self.SwitchAnalysis.process_ensemble_variation, args=[new_shuttle])
            mew_thread.start()

    
    def wait_for_finish(self):
                
        stop:bool = False
        while stop == False:
            print(f'Checking LMV status at {datetime.now()}')
            current_status: List[bool] = self.group_token.group_done_status
            is_done:bool = self.group_token.is_done
            
            message: str = ''
            for index in range(self.num_groups):
                goup_value:str = self.group_token.group_values[index]
                done_status: bool = self.group_token.group_done_status[index]
                message = message + f'Group_{index+1}: kcal_group={goup_value}, status={done_status}\n'
            print(message)

            if is_done == True:
                stop = True
                print(f'Its done at {datetime.now()}')
            else:
                dwell_time:int = 5
                print(f'dwelling for {dwell_time} seconds until next check')
                time.sleep(dwell_time)


class LMV():

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

    class LMV_ThreadProcessor():
        
        def __init__(self, stuctures: List[Sara2StructureList], comp_structure: Sara2SecondaryStructure = Sara2SecondaryStructure(), comp_struct_list_option: List[Sara2SecondaryStructure] = []) -> None:
            self._sara2_groups: List[Sara2StructureList] = stuctures
            num_groups:int = len(stuctures)
            self._num_groups: int =  num_groups
            self._group_token: LMV_Token = LMV_Token(num_groups)
            self._LMV: EnsembleVariation = EnsembleVariation()
            self._comparison_structure: Sara2SecondaryStructure = comp_structure

        @property
        def sara2_groups(self):
            return self._sara2_groups

        @sara2_groups.setter
        def sara2_groups(self, new_list:List[Sara2StructureList]):
            self._sara2_groups = new_list
        
        @property
        def comparison_structure(self):
            return self._comparison_structure

        @comparison_structure.setter
        def comparison_structure(self, new_struct:Sara2SecondaryStructure):
            self._comparison_structure = new_struct
        
        @property
        def comp_struct_list_option(self):
            return self._comp_struct_list_option

        @comp_struct_list_option.setter
        def comp_struct_list_option(self, new_list:List[Sara2SecondaryStructure]):
            self._comp_struct_list_option = new_list

        @property
        def num_groups(self):
            return self._num_groups

        @num_groups.setter
        def num_groups(self, new_num:int):
            self._num_groups = new_num

        @property
        def group_token(self):
            return self._group_token

        @group_token.setter
        def group_token(self, new_token:LMV_Token):
            self._group_token = new_token
        
        @property
        def LMV(self):
            return self._LMV

        @LMV.setter
        def LMV(self, new_lmv:EnsembleVariation):
            self._LMV = new_lmv

        def run_LMV(self):
            self.start_calculations()
            self.wait_for_finish()
            #the test should be done now
            #check foor index that is -1 and if so then use prev value
            num_groups:int = len(self.group_token.group_results)
            for index in range(1, num_groups):
                if self.group_token.group_results[index].ev_normalized == -1:
                    previous_EV = self.group_token.group_results[index-1]
                    self.group_token.group_results[index] = previous_EV
                    self.group_token.group_dict[index] = previous_EV
            return self.group_token

        def start_calculations(self):
            comp_structure: Sara2SecondaryStructure = Sara2SecondaryStructure()               
            for thread_index in range(self.num_groups):
                if len(self.comp_struct_list_option) == self.num_groups:
                    comp_structure = self.comp_struct_list_option[thread_index]
                else:
                    comp_structure = self.comparison_structure
                sara2_structs: Sara2StructureList  = self.sara2_groups[thread_index]
                new_shuttle: LMV_Shuttle = LMV_Shuttle(structs_list=sara2_structs, mfe=comp_structure, group_index=thread_index,token=self.group_token) 
                mew_thread = threading.Thread(target=self.LMV.thread_EV, args=[new_shuttle])
                mew_thread.start()

        
        def wait_for_finish(self):
                    
            stop:bool = False
            while stop == False:
                print(f'Checking LMV status at {datetime.now()}')
                current_status: List[bool] = self.group_token.group_done_status
                is_done:bool = self.group_token.is_done
                
                message: str = ''
                for index in range(self.num_groups):
                    goup_value:str = self.group_token.group_values[index]
                    done_status: bool = self.group_token.group_done_status[index]
                    message = message + f'Group_{index+1}: kcal_group={goup_value}, status={done_status}\n'
                print(message)

                if is_done == True:
                    stop = True
                    print(f'Its done at {datetime.now()}')
                else:
                    dwell_time:int = 5
                    print(f'dwelling for {dwell_time} seconds until next check')
                    time.sleep(dwell_time)