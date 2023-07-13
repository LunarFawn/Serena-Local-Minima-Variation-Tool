
from dataclasses import dataclass
from typing import List, Dict
from enum import Enum
from threading import Thread
import threading
from time import sleep

from serena.apps.original_weighted_analysis import Sara2SecondaryStructure, Sara2StructureList, EnsembleVariation, EVResult
from serena.utilities.results import SwitchAnalysisResults
from serena.utilities.data import SwitchAnalysisDataToProcess



class CPU_Data():
    
    def __init__(self) -> None:
        self._temperature_data:SwitchAnalysisDataToProcess = SwitchAnalysisDataToProcess()
    
    @property
    def temperature_data(self):
        return self._temperature_data

    @temperature_data.setter
    def temperature_data(self, temperature_data:SwitchAnalysisDataToProcess):
        self.temperature_data = temperature_data
    

class CPU_Result():    
    def __init__(self) -> None:
        self._temperature_result:SwitchAnalysisResults = SwitchAnalysisResults()

    @property
    def temperature_result(self):
        return self._temperature_result

    @temperature_result.setter
    def temperature_result(self, temperature_result:SwitchAnalysisResults):
        self._temperature_result = temperature_result
    
  
class CPU_Token():
    """
    This is the class that contains the token that
    holds all the results of the functions on the threads
    """
    def __init__(self, total_threads: int) -> None:
        self._thread_results: List[CPU_Result] = total_threads * [CPU_Result()]
        self._thread_dict: Dict[int,CPU_Result] = {}
        self._thread_comments: List[str] = total_threads * ['']
        self._thread_done_status: List[bool] = total_threads * [False]
    
    @property
    def thread_results(self):
        return self._thread_results
        
    def set_thread_results(self, index:int, value:CPU_Result):
        self._thread_results[index]=value

    @property
    def thread_dict(self):
        return self._thread_dict
        
    def set_thread_dict(self, index:int, value:CPU_Result):
        self._thread_dict[index]=value
    
    @property
    def thread_comments(self):
        return self._thread_comments
        
    def set_thread_comments(self, index:int, value:str):
        self._thread_comments[index]=value

    @property
    def thread_done_status(self):
        return self._thread_done_status
        
    def set_thread_done_status(self, index:int, state:bool):
        self._thread_done_status[index]=state
    
    @property
    def is_done(self):
        is_completed:bool = False
        if self._thread_done_status.count(False) == 0:
            #its done
            is_completed = True
        return is_completed


class LMV_Shuttle():
    """
    This is the class the transports the data to be processed to thefunction in the thread
    from the original caller of the threads
    """
    def __init__(self, data_to_process:CPU_Data, thread_index:int, token:CPU_Token) -> None:
        self._data_to_process:CPU_Data = data_to_process
        self._thread_index:int = thread_index
        self._token:CPU_Token = token
    
    @property
    def data_to_process(self):
        return self._data_to_process

    @data_to_process.setter
    def data_to_process(self, data_to_process: CPU_Data):
        self._data_to_process = data_to_process
    
    @property
    def thread_index(self):
        return self._thread_index

    @thread_index.setter
    def thread_index(self, thread_index: int):
        self._thread_index = thread_index
    
    @property
    def token(self):
        return self._token

    @token.setter
    def token(self, new_token: CPU_Token):
        self._token = new_token


class CPU_ThreadProcessor():
    """
    This is the class that actuallt runs the threads and is called from the main caller
    """  
    def __init__(self, data_to_process:CPU_Data, function_to_run, total_threads:int) -> None:
        self._data_to_process:CPU_Data = data_to_process
     
        self._total_threads: int =  total_threads
        self._thread_token: CPU_Token = CPU_Token(total_threads)
        self._function_to_run = function_to_run
    
    @property
    def data_to_process(self):
        return self._data_to_process

    @data_to_process.setter
    def data_to_process(self, data_to_process:CPU_Data):
        self._data_to_process = data_to_process

    @property
    def total_threads(self):
        return self._total_threads

    @total_threads.setter
    def total_threads(self, total_threads:int):
        self._total_threads = total_threads
    
    @property
    def thread_token(self):
        return self._thread_token

    @thread_token.setter
    def thread_token(self, thread_token:CPU_Token):
        self._thread_token = thread_token

    @property
    def function_to_run(self):
        return self._function_to_run

    @function_to_run.setter
    def function_to_run(self, function_to_run):
        self._function_to_run = function_to_run
    
    def run_cpu(self):
        self.start_functions()
        self.wait_for_finish()
        return self.thread_token

    def start_functions(self):
        for thread_index in range(self.total_threads):
            new_shuttle: LMV_Shuttle = LMV_Shuttle(data_to_process=self.data_to_process, thread_index=thread_index, token=self.thread_token) 
            mew_thread = Thread(target=self.function_to_run, args=[new_shuttle])
            mew_thread.start()

    
    def wait_for_finish(self):
                
        stop:bool = False
        while stop == False:
            #print(f'Checking LMV status at {datetime.now()}')
            current_status: List[bool] = self.thread_token.thread_done_status
            is_done:bool = self.thread_token.is_done
            
            message: str = ''
            for index in range(self.total_threads):
                thread_comment:str = self.thread_token.set_thread_comments[index]
                done_status: bool = self.thread_token.thread_done_status[index]
                message = message + f'Thread_{index+1}, status={done_status}, comment={thread_comment}\n'
            #print(message)

            if is_done == True:
                stop = True
                #print(f'Its done at {datetime.now()}')
            else:
                dwell_time:int = 5
                #print(f'dwelling for {dwell_time} seconds until next check')
                sleep(dwell_time)
        



