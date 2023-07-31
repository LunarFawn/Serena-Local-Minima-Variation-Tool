"""
File to hold the main entry point for determining if a rna
sequence's ensemble is accepting of a switch
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
from serena.utilities.ensemble_analysis import InvestigateEnsemble
from serena.utilities.nupack4_sara2_extension import NUPACK4Interface, MaterialParameter, NupackSettings
from serena.utilities.weighted_structures import WeightedStructure




class SwitchAccetance():

    def __init__(self) -> None:
        pass

    def rna_switch_scoring(self):
        material:MaterialParameter = MaterialParameter.rna95_nupack4
        temp:int = 37
        kcal_span:int = 7
        kcal_unit_increments:float = 1
        sequence:str = ''

    def run_eterna_pnas(self):

        vienna2_fmn_hack: Vienna2FMNInterface = Vienna2FMNInterface()

        details:str= 'unbound_as_target_run4'#f'20k_filtered_weighted_100K_gtrequal2_nucpenalty_run_1ish'
        pnas_round101_sheet:str = 'R101 Filtered Switch All'
        same_state:str='1'
        sublab_name:str = f'Same State NG {same_state}'
        run_name:str = f'SSNG{same_state}_{details}'


        pnas_path:str = '/mnt/g/serena/pnas.2112979119.sd01_eternacon.xlsx'
        timestr = time.strftime("%Y%m%d-%H%M%S")
        save_path:str = f'/mnt/g/serena/{run_name}/pnas.2112979119.sd01_eternacon_{timestr}.xlsx'
        
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
                fmn_weighted_struct: WeightedStructure = WeightedStructure()
                struct_to_use= fmn_weighted_struct.make_weighted_struct(fmn_subopt)
            else:
                fmn_struct = vienna2_fmn_hack.rnafold_fmn(input_sequence=sequence,
                                                            do_fmn=True)
                struct_to_use = fmn_struct.structure

            
