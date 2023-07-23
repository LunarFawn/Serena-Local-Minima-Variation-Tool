"""
Stuff I spun up just for eternacon presentation
"""


from dataclasses import dataclass
import os
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
from serena.utilities.logging_this import PNASAnalysisLogging

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

        details:str= 'unbound_as_target_run4'#f'20k_filtered_weighted_100K_gtrequal2_nucpenalty_run_1ish'
        pnas_round101_sheet:str = 'R101 Filtered Switch All'
        same_state:str='1'
        sublab_name:str = f'Same State NG {same_state}'
        run_name:str = f'SSNG{same_state}_{details}'


        switch:OriginalSwitchAnalysis = OriginalSwitchAnalysis()
        switch.save_folder_path = f'/mnt/g/serena/{run_name}'
        switch.sublab_name = f'SSNG{same_state}'
        

        pnas_path:str = '/mnt/g/serena/pnas.2112979119.sd01_eternacon.xlsx'
        timestr = time.strftime("%Y%m%d-%H%M%S")
        save_path:str = f'/mnt/g/serena/{run_name}/pnas.2112979119.sd01_eternacon_{timestr}.xlsx'
        
        
        

        if os.path.exists(switch.save_folder_path) is False:
            os.mkdir(switch.save_folder_path)
            time.sleep(1)
        
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

            #make a new line of just this designs row
            
            
            do_weighted:bool = False
            struct_to_use:str = ''
            if do_weighted is True:
            #this is the fmn bound mfe struct, subopt list and weighted struck
                fmn_subopt = vienna2_fmn_hack.rnasubopt_fmn(input_sequence=sequence,
                                                            do_fmn=True)
                fmn_weighted_struct: WeightedStructure = WeightedStructure(fmn_subopt)
                struct_to_use= fmn_weighted_struct.weighted_structure.structure
            else:
                fmn_struct = vienna2_fmn_hack.rnafold_fmn(input_sequence=sequence,
                                                            do_fmn=True)
                struct_to_use = fmn_struct.structure
            

            analysis:PredictionReponse = switch.do_switch_analysis(sequence=sequence,
                                                        fmn_struct=struct_to_use,
                                                        fmn_struct_free_energy=0,
                                                        span=7,
                                                        units=1,
                                                        run_name=design_id,
                                                        manual=False)

            pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'PredictedFoldChange'] = analysis.foldchange
            #design_data_df['PredictedFoldChange'] = analysis.foldchange

            pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'AvgSwitchScore'] = statistics.fmean(analysis.raw_scores)
            #design_data_df['AvgSwitchScore'] = statistics.fmean(analysis.raw_scores)
            pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'AvgNumSstruct'] = statistics.fmean(analysis.num_structs)
            #design_data_df['AvgNumSstruct'] = statistics.fmean(analysis.num_structs)
            
            pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, '36Deg_SwitchScore'] = analysis.raw_scores[0]
            #design_data_df['36Deg_SwitchScore'] = analysis.raw_scores[0]
            pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, '37Deg_SwitchScore'] = analysis.raw_scores[1]
            #design_data_df['37Deg_SwitchScore'] = analysis.raw_scores[1]
            pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, '38Deg_SwitchScore'] = analysis.raw_scores[2]
            #design_data_df['38Deg_SwitchScore'] = analysis.raw_scores[2]

            pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, '36Deg_NumStructs'] =  analysis.num_structs[0]
            #design_data_df['36Deg_NumStructs'] = analysis.num_structs[0]
            pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, '37Deg_NumStructs'] =  analysis.num_structs[1]
            #design_data_df['37Deg_NumStructs'] = analysis.num_structs[1]
            pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, '38Deg_NumStructs'] =  analysis.num_structs[2]
            #design_data_df['38Deg_NumStructs'] = analysis.num_structs[2]
            
            design_data_df:DataFrame = pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID]
            logging: PNASAnalysisLogging = PNASAnalysisLogging()
            logging.save_excel_sheet(design_data_df, save_path, sublab_name)
           
     
        print("Its done!!!")

eternacon_stuff:Eternacon2023 = Eternacon2023()
eternacon_stuff.run_nupack_vienna()

