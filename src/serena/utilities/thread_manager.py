
from dataclasses import dataclass
from typing import List, Dict
from enum import Enum

from serena.apps.original_weighted_analysis import Sara2SecondaryStructure, Sara2StructureList, EnsembleVariation, EVResult





@dataclass
class SwitchInput():
    sequence:str
    fmn_struct:str
    fmn_struct_free_energy:float
    span:int
    units:int
    run_name:str

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

