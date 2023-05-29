"""
Class for generating a report file for each analysis type
"""
from pathlib import Path
import time
from typing import List, Dict
from dataclasses import dataclass

import serena.ensemble_variation as ev
from serena.ensemble_variation import EnsembleVariation, EVResult
import serena.structures as ser_structs
from serena.structures import MultipleEnsembleGroups
import serena.nupack4_sara2_extension as nupack_extension
from serena.nupack4_sara2_extension import NUPACK4Interface, NupackSettings, MaterialParameter


class LocalMinimaVariationReport():
    """
    Class for making the local minima variation (ensemble variation) report
    """

    class LocalMinimaVariationData():
        def __init__(self) -> None:
            pass        


    def __init__(self, working_folder: Path) -> None:
        self._working_folder: Path = working_folder
    
    @property
    def working_folder(self):
        return self._working_folder
    
    @working_folder.setter
    def working_folder(self, folder_path:Path):
        self._working_folder = folder_path

    def generate_text_report(self, report_name:str):
        #now save data to csv file        
        timestr = time.strftime("%Y%m%d-%H%M%S")

        time_span: List[float] = []
        tick_span = []
        index_span = range(len(ev_result_mfe.groups_list))
        

        mfe_value:float = ev_result_mfe.groups_list[0].mfe_freeEnergy
        seed_value:float = mfe_value
        tick_value:float = 0
        
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

        new_switch_string: List[float] = []
        for ev in switch_result_target.group_ev_list:
            ev_value = ev.ev_normalized
            new_switch_string.append(ev_value)

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
            LMSV_US_target = new_switch_string[index]
            LMSV_US_folded = new_switch_string_folded[index]
            line:str = f'{kcal},{LMSV_U_mfe},{LMSV_U_rel},{LMSV_US_target},{LMSV_US_folded}\n'
            csv_log_results.append(line)

        csv_record_pathstr = f'{self._working_folder}/{report_name}_{timestr}.csv'
        csv_lines:List[str]=[]
        with open(csv_record_pathstr, 'w') as csv_file:
            #first write teh header
            csv_lines.append(f'Local Minima Structure Variation Data\n')
            csv_lines.append(f'Creation Date={datetime.now()}\n')
            csv_lines.append("---------------------------------------\n")
            csv_lines.append("***DESIGN INFO***\n")
            csv_lines.append(f'Design Name = {name}\n')
            csv_lines.append(f'DesignID = {designID}\n')
            csv_lines.append(f'Lab Name = {labname}\n')
            csv_lines.append(f'Sequence = {sequence}\n')
            csv_lines.append(f'Eterna_Score = {eterna_score}\n')
            csv_lines.append(f'FoldChange = {fold_change}\n')
            csv_lines.append(f'2nd State Target Structure = {target}\n')
            csv_lines.append(f'2nd State Folded Structure = {folded}\n')
            csv_lines.append(f'2nd State Folded Oligo Energy = {folded_energy_ligoligo}\n')
            csv_lines.append(f'Energy Span from MFE = {span}\n')
            csv_lines.append(f'Energy span units = {units}\n')
            csv_lines.append("---------------------------------------\n")
            csv_lines.append("***RAW DATA***\n")
            csv_lines = csv_lines + csv_log_results
            csv_lines.append("---------------------------------------\n")
            csv_lines.append("***METRICS***\n")
            csv_lines.append(f'Polymorphicity Level (2kcal to end of sample) = {mfe_fold_EV_delta}\n')
            csv_lines.append(f'LMV_US_folded at folded with lignad/oligo energy = {ev_oligo_folded}\n')
            csv_lines.append("---------------------------------------\n")
            csv_lines.append("EOF\n")
            csv_file.writelines(csv_lines)
    

