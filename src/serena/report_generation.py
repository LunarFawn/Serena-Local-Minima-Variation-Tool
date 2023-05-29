"""
Class for generating a report file for each analysis type
"""
from pathlib import Path
import time
from typing import List, Dict
from dataclasses import dataclass
from datetime import datetime

import serena.ensemble_variation as ev
from serena.ensemble_variation import EnsembleVariation, EVResult
import serena.structures as ser_structs
from serena.structures import MultipleEnsembleGroups
import serena.nupack4_sara2_extension as nupack_extension
from serena.nupack4_sara2_extension import NUPACK4Interface, NupackSettings, MaterialParameter

@dataclass
class SequenceResults():
    result_data: EVResult
    result_name: str
    
@dataclass
class SequenceInfo():
    lab_name:str = ''
    sequence: str = ''
    sequence_name: str = ''
    sequence_ID: int = -1
    folded_energy: float = 0
    ligand_oligo_name:str = ''
    eterna_score:float = 0
    fold_change:float = 0
    number_of_clusters:int = 0
    temp_C: int = 0
    span = 0
    units = 0
    folded_struct:str


class FullRunInfo():

        def __init__(self, run_name:str, run_ID: int, sequence_info: SequenceInfo) -> None:
            self._run_results: List[SequenceResults] = []
            self._sequence_info: SequenceInfo  = sequence_info
            self._run_name: str = run_name
            self._run_ID: int = run_ID
        
        def add_sequence_result(self, result: SequenceResults):
            self._run_results.append(result)
        
        @property
        def run_results(self):
            return self._run_results
        
        @run_results.setter
        def run_results(self, results:List[SequenceResults] ):
            self._run_results = results

        @property
        def sequence_info(self):
            return self._sequence_info
        
        @sequence_info.setter
        def sequence_info(self, info:SequenceInfo):
            self._sequence_info = info
        
        @property
        def run_name(self):
            return self._run_name
        
        @run_name.setter
        def run_name(self, name:str):
            self._run_name = name

        @property
        def run_ID(self):
            return self._run_ID
        
        @run_ID.setter
        def run_ID(self, ID:int):
            self._run_ID = ID


class LocalMinimaVariationReport():
    """
    Class for making the local minima variation (ensemble variation) report
    """
    
    def __init__(self, working_folder: Path) -> None:
        self._working_folder: Path = working_folder
    
    @property
    def working_folder(self):
        return self._working_folder
    
    @working_folder.setter
    def working_folder(self, folder_path:Path):
        self._working_folder = folder_path

    def generate_text_report(self, run_info: FullRunInfo):
        #now save data to csv file        
        timestr = time.strftime("%Y%m%d-%H%M%S")

        ev_result_mfe: EVResult = run_info.run_results[0].result_data
        ev_result_rel: EVResult = run_info.run_results[1].result_data
        switch_result_folded: EVResult = run_info.run_results[2].result_data

        time_span: List[float] = []
        tick_span = []
        index_span = range(len(ev_result_mfe.groups_list))
        

        mfe_value:float = ev_result_mfe.groups_list[0].mfe_freeEnergy
        seed_value:float = mfe_value
        tick_value:float = 0
        
        units:int = run_info.sequence_info.units

        num_samples: int = len(ev_result_mfe.groups_list)
        for index in range(num_samples):
            seed_value = seed_value + float(units)
            tick_value = tick_value + float(units)
            #time_span is teh MFE values
            #tick_span is the index value (i.e. 0, 0.5, 1, 1.5)
            time_span.append(seed_value)
            tick_span.append(tick_value)

        

        new_list_string_mfe: List[float] = []
        for ev in ev_result_mfe.group_ev_list:
            ev_value = ev.ev_normalized
            new_list_string_mfe.append(ev_value)

        new_list_string_rel: List[float] = []
        for ev in ev_result_rel.group_ev_list:
            ev_value = ev.ev_normalized
            new_list_string_rel.append(ev_value)

        new_switch_string_folded: List[float] = []
        for ev in switch_result_folded.group_ev_list:
            ev_value = ev.ev_normalized
            new_switch_string_folded.append(ev_value)

        csv_log_results: List[str]=[]
        csv_log_results.append("Kcal,LMSV_U_mfe,LMSV_U_rel,LMSV_US_target,LMSV_US_folded\n")
        for index in range(len(new_list_string_mfe)):
            kcal = time_span[index]
            LMSV_U_mfe = new_list_string_mfe[index]
            LMSV_U_rel = new_list_string_rel[index]
            LMSV_US_folded = new_switch_string_folded[index]
            line:str = f'{kcal},{LMSV_U_mfe},{LMSV_U_rel},{LMSV_US_folded}\n'
            csv_log_results.append(line)


        csv_record_pathstr = f'{self._working_folder}/{run_info.run_name}_{timestr}.csv'
        csv_lines:List[str]=[]
        with open(csv_record_pathstr, 'w') as csv_file:
            #first write teh header
            csv_lines.append(f'Local Minima Structure Variation Data\n')
            csv_lines.append(f'Creation Date={datetime.now()}\n')
            csv_lines.append("---------------------------------------\n")
            csv_lines.append("***DESIGN INFO***\n")
            csv_lines.append(f'Design Name = {run_info.sequence_info.sequence_name}\n')
            csv_lines.append(f'DesignID = {run_info.sequence_info.sequence_ID}\n')
            csv_lines.append(f'Lab Name = {run_info.sequence_info.lab_name}\n')
            csv_lines.append(f'Sequence = {run_info.sequence_info.sequence}\n')
            csv_lines.append(f'Eterna_Score = {run_info.sequence_info.eterna_score}\n')
            csv_lines.append(f'FoldChange = {run_info.sequence_info.fold_change}\n')
            csv_lines.append(f'2nd State Folded Structure = {run_info.sequence_info.folded_struct}\n')
            csv_lines.append(f'2nd State Folded Oligo Energy = {run_info.sequence_info.folded_energy}\n')
            csv_lines.append(f'Energy Span from MFE = {run_info.sequence_info.span}\n')
            csv_lines.append(f'Energy span units = {units}\n')
            csv_lines.append("---------------------------------------\n")
            csv_lines.append("***RAW DATA***\n")
            csv_lines = csv_lines + csv_log_results
            csv_lines.append("---------------------------------------\n")
            csv_lines.append("EOF\n")
            csv_file.writelines(csv_lines)
    

