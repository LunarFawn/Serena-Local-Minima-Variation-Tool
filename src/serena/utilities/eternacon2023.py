"""
Stuff I spun up just for eternacon presentation
"""


from dataclasses import dataclass
import statistics
import openpyxl
import pandas as pd
from pandas import DataFrame
from typing import List, Dict
import time


from serena.utilities.Sara2_API_Python3 import Sara2API, puzzleData
from serena.utilities.vienna2_fmn_hack_interface import Vienna2FMNInterface
from serena.utilities.ensemble_structures import Sara2SecondaryStructure, Sara2StructureList, ComparisonStructures
from serena.utilities.weighted_structures import WeightedStructure
from serena.scripts.run_switch_analysis import OriginalSwitchAnalysis, PredictionReponse
from serena.utilities.logging import PNASAnalysisLogging

@dataclass
class PNASWetLab():
    pass

class Eternacon2023():

    def __init__(self) -> None:
        pass

    def run_nupack_vienna(self):
        """
        run this
        """
        
        vienna2_fmn_hack: Vienna2FMNInterface = Vienna2FMNInterface()

        pnas_path:str = '/mnt/g/serena/pnas.2112979119.sd01_eternacon.xlsx'
        timestr = time.strftime("%Y%m%d-%H%M%S")
        save_path:str = f'/mnt/g/serena/test/pnas.2112979119.sd01_eternacon_{timestr}.xlsx'
        pnas_round101_sheet:str = 'Round 7 (R101) (2)'
        sublab_name:str = 'Same State NG 1'
        
        new_sara:Sara2API = Sara2API()
        puzzle_data: puzzleData
        pandas_sheet: DataFrame
        puzzle_data, pandas_sheet = new_sara.ProcessLab(path=pnas_path,
                                                      designRound_sheet=pnas_round101_sheet,
                                                      sublab_name=sublab_name)


        predicted_foldchange_list: List[float] = []
        avg_raw_score_list: List[float] = []
        avg_num_structures_list: List[float] = []

        raw_score_36_list: List[float] = []
        raw_score_37_list: List[float] = []
        raw_score_38_list: List[float] = []

        num_structures_36_list: List[float] = []
        num_structures_37_list: List[float] = []
        num_structures_38_list: List[float] = []

        flag:int =0
        for design in puzzle_data.designsList:
            design_id= str(design.design_info.DesignID)
            sequence = design.design_info.Sequence
            fold_change = design.wetlab_results.FoldChange
            eterna_score = design.wetlab_results.Eterna_Score
            folding_subscore = design.wetlab_results.Folding_Subscore
            switch_subscore = design.wetlab_results.Switch_Subscore
            baseline_subscore = design.wetlab_results.Baseline_Subscore

            #this is the fmn bound mfe struct, subopt list and weighted struck
            fmn_struct = vienna2_fmn_hack.rnafold_fmn(sequence)
            fmn_subopt = vienna2_fmn_hack.rnasubopt_fmn(sequence)
            fmn_weighted_struct: WeightedStructure = WeightedStructure(fmn_subopt)

            switch:OriginalSwitchAnalysis = OriginalSwitchAnalysis()
            switch.save_folder_path = '/mnt/g/serena/test'
            switch.sublab_name = "SSNG1"

            analysis:PredictionReponse = switch.do_switch_analysis(sequence=sequence,
                                                        fmn_struct=fmn_struct.structure,
                                                        fmn_struct_free_energy=0,
                                                        span=7,
                                                        units=1,
                                                        run_name=design_id,
                                                        manual=False)
            
            predicted_foldchange_list.append(analysis.foldchange)
            
            
            avg_raw_score_list.append(statistics.fmean(analysis.raw_scores))
            avg_num_structures_list.append(statistics.fmean(analysis.num_structs))
            
            raw_score_36_list.append(design_id)
            raw_score_37_list.append(design_id)
            raw_score_38_list.append(design_id)
            
            num_structures_36_list.append(design_id)
            num_structures_37_list.append(design_id)
            num_structures_38_list.append(design_id)
            
        pandas_sheet['PredictedFoldChange'] = predicted_foldchange_list
        
        pandas_sheet['AvgSwitchScore'] = avg_raw_score_list
        pandas_sheet['AvgNumSstruct'] = avg_num_structures_list
        
        pandas_sheet['36Deg_SwitchScore'] = raw_score_36_list
        pandas_sheet['37Deg_SwitchScore'] = raw_score_37_list
        pandas_sheet['38Deg_SwitchScore'] = raw_score_38_list

        pandas_sheet['36Deg_NumStructs'] = num_structures_36_list
        pandas_sheet['37Deg_NumStructs'] = num_structures_37_list
        pandas_sheet['38Deg_NumStructs'] = num_structures_38_list

        logging: PNASAnalysisLogging = PNASAnalysisLogging()
        logging.save_dataframe_to_excel(src_dataframe=pandas_sheet,
                                        dst_path=save_path,
                                        dst_sheet_name=sublab_name)
     
        print("Its done!!!")

eternacon_stuff:Eternacon2023 = Eternacon2023()
eternacon_stuff.run_nupack_vienna()

