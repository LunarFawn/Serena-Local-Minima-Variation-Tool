"""
File for handeling generation of reports on the investigator results as well as 
handeling the console script entry points for it
"""
from pathlib import Path
import time
from typing import List, Dict
from dataclasses import dataclass
from datetime import datetime
import numpy as np
from pandas import DataFrame
import os
from matplotlib.markers import MarkerStyle
import argparse

from enum import Enum

import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import matplotlib

from serena.utilities.ensemble_variation import EnsembleVariation, EVResult

from serena.analysis.investigator import (
    InvestigatorResults,
    RatioResults,
    
    
)

from serena.bin.backup_serena_v2 import ArchiveSecondaryStructureList

from serena.utilities.comparison_structures import (
    ComparisonNucCounts,
    ComparisonNucResults,
)

from serena.utilities.ensemble_structures import (
    Sara2SecondaryStructure,Sara2StructureList
)

from serena.utilities.ensemble_variation import (
    EV
)

from serena.utilities.local_minima_variation import (
    ComparisonLMV,
    ComparisonLMVResponse,
    
)

from serena.interfaces.Sara2_API_Python3 import (
    DesignPerformanceData,
    DesignInformation,
    WetlabData,
    Sara2API,
    puzzleData
)

from serena.analysis.ensemble_analysis import InvestigateEnsembleResults

from serena.bin.backup_investigator_v1 import ArchiveInvestigator
from serena.scripts.analyze_pnas_2112979119_sd01 import ProcessPNAS, ArchiveInvestigatorData, ArchiveFlow

from serena.utilities.structure_search import StaticSystemDetector, PrimeNucCounts, MolecularSnareDetector, MoleculareSnareDef, SnareResults

from serena.scripts.analyze_pnas_2112979119_sd01 import ArchiveData, ArchiveFlow

class archiveType(Enum):
    RATIO='RATIO'
    COUNT="COUNT"
    LMV="LMV"
    LMV_REL="LMV_REL"
    LMV_MFE="LMV_MFE"
    LMV_COMP="LMV_COMP"
    STATIC_PRIMES="STATIC_PRIMES"
    SNARE='SNARE'
    
class ScoreType(Enum):
    BASELINE='BASELINE'
    FOLDING='FOLDING'
    SWITCH='SWITCH'
    FOLDCHANGE='FOLDCHANGE'
    KDON='KDON'
    KDOFF='KDOFF'
    ETERNA="ETERNA"
    BASIC="BASIC"
    ADVANCED="ADVANCED"
    
    

class InvestigatorReportGeneration():
    
    def __init__(self) -> None:
        pass
    
            
    def generate_nuc_count_plot(self, timestr:str, data:List[ArchiveInvestigatorData], source_data:List[ArchiveData], attr:archiveType, nuc_count_name:str, x_string:str, x_range:float, training:bool = True, score_type:ScoreType=ScoreType.FOLDCHANGE, snare_binding:str=None ):
        
        
        #find how many enesmble energy groups there are
        num_groups:int = data[0].investigator.investigator_results.num_groups
        
        
        
        ax:plt = None
        fig, ax = plt.subplots(num_groups, constrained_layout=True, figsize=(15, 15))
        fig.suptitle(nuc_count_name)
        fig.supxlabel(x_string)
        
        
        # fig.text(0.50, 0.02, 
        #      "Count of nucleotide with same pairing in unbound and bound state", 
        #     horizontalalignment='center', wrap=True )
        
        for plt_index in range(num_groups):
            
            #get all the data for each plot
            
            #first get the both count
            good_list:List[float] = []
            good_fold_change:List[float] = []
            bad_list:List[float] = []
            bad_change:List[float] = []
            
            result_list:List[float] = []
            
            low_structs_num:int = 4000#5000#3000
            med_structs_num:int = 8000#15000#6000
            high_structs_num:int = 12000#25000#9000
            
            low_structs:List[float] = []
            med_structs:List[float] = []
            high_structs:List[float] = []
            obsurd_structs:List[float] = []
            
            
            baseline_scores:List[float] = []
            baseline_low_structs:List[float] = []
            baseline_med_structs:List[float] = []
            baseline_high_structs:List[float] = []
            baseline_obsurd_structs:List[float] = []
                        
            folding_scores:List[float] = []
            folding_low_structs:List[float] = []
            folding_med_structs:List[float] = []
            folding_high_structs:List[float] = []
            folding_obsurd_structs:List[float] = []
            
            switch_scores:List[float] = []
            switch_low_structs:List[float] = []
            switch_med_structs:List[float] = []
            switch_high_structs:List[float] = []
            switch_obsurd_structs:List[float] = []
            
            kdoff_values:List[float] = []
            kdoff_low_structs:List[float] = []
            kdoff_med_structs:List[float] = []
            kdoff_high_structs:List[float] = []
            kdoff_obsurd_structs:List[float] = []
            
            kdone_values:List[float] = []
            kdone_low_structs:List[float] = []
            kdone_med_structs:List[float] = []
            kdone_high_structs:List[float] = []
            kdone_obsurd_structs:List[float] = []
            
            eterna_values:List[float] = []
            eterna_low_structs:List[float] = []
            eterna_med_structs:List[float] = []
            eterna_high_structs:List[float] = []
            eterna_obsurd_structs:List[float] = []
            
            serena_basic_values:List[float] = []
            serena_basic_low_structs:List[float] = []
            serena_basic_med_structs:List[float] = []
            serena_basic_high_structs:List[float] = []
            serena_basic_obsurd_structs:List[float] = []
            
            serena_advanced_values:List[float] = []
            serena_advanced_low_structs:List[float] = []
            serena_advanced_med_structs:List[float] = []
            serena_advanced_high_structs:List[float] = []
            serena_advanced_obsurd_structs:List[float] = []
            
            result_fold_change:List[float] = []
            result_fold_change_low_structs:List[float] = []
            result_fold_change_med_structs:List[float] = []
            result_fold_change_high_structs:List[float] = []
            result_fold_change_obsurd_structs:List[float] = []
            
            for index, design in enumerate(data):
                if training is True:
                    if design.design_info.design_info.Puzzle_Name == 'good':
                        if attr == archiveType.COUNT:
                            good_list.append(getattr(design.investigator.investigator_results.comp_nuc_counts.comparison_nuc_counts[plt_index], nuc_count_name)) 
                        elif attr == archiveType.RATIO:
                            good_list.append(getattr(design.investigator.investigator_results.comparison_eval_results.ratios[plt_index], nuc_count_name)) 
                        elif attr == archiveType.LMV:
                            new_attr:EV = getattr(design.investigator.investigator_results.lmv_values.lmv_comps[plt_index], nuc_count_name)
                            good_list.append(new_attr.ev_normalized)
                        # elif attr == archiveType.LMV_MFE:
                        #     good_list.append(getattr(design.investigator.investigator_results.lmv_values.lmv_comps[plt_index], nuc_count_name).ev_normalized)
                        # elif attr == archiveType.LMV_COMP:
                        #     good_list.append(getattr(design.investigator.investigator_results.lmv_values.lmv_comps[plt_index], nuc_count_name).ev_normalized)    
                        good_fold_change.append(design.design_info.wetlab_results.FoldChange)
                    elif design.design_info.design_info.Puzzle_Name == 'bad':
                        if attr == archiveType.COUNT:
                            bad_list.append(getattr(design.investigator.investigator_results.comp_nuc_counts.comparison_nuc_counts[plt_index], nuc_count_name)) 
                        elif attr == archiveType.RATIO:
                            bad_list.append(getattr(design.investigator.investigator_results.comparison_eval_results.ratios[plt_index], nuc_count_name))  
                        elif attr == archiveType.LMV:
                            new_attr:EV = getattr(design.investigator.investigator_results.lmv_values.lmv_comps[plt_index], nuc_count_name)
                            bad_list.append(new_attr.ev_normalized)
                        # elif attr == archiveType.LMV_MFE:
                        #     bad_list.append(getattr(design.investigator.investigator_results.lmv_values.lmv_comps[plt_index], nuc_count_name).ev_normalized)
                        # elif attr == archiveType.LMV_COMP:
                        #     bad_list.append(getattr(design.investigator.investigator_results.lmv_values.lmv_comps[plt_index], nuc_count_name).ev_normalized)  
                        bad_change.append(design.design_info.wetlab_results.FoldChange)
                else:
                    new_attr_value:float
                    if attr == archiveType.COUNT:
                        new_attr_value = getattr(design.investigator.investigator_results.comp_nuc_counts.comparison_nuc_counts[plt_index], nuc_count_name)
                        # if design.investigator.number_structures <= low_structs:
                        #     low_structs.append(count_attr)
                        # elif design.investigator.number_structures > low_structs_num and design.investigator.number_structures <= med_structs_num:
                        #     med_structs.append(count_attr)
                        # elif design.investigator.number_structures > med_structs_num and design.investigator.number_structures <= high_structs_num:
                        #     high_structs.append(count_attr)
                        # elif design.investigator.number_structures > high_structs_num:
                        #     obsurd_structs.append(count_attr) 
                        
                            
                        # result_list.append(count_attr) 
                    elif attr == archiveType.RATIO:
                        new_attr_value = getattr(design.investigator.investigator_results.comparison_eval_results.ratios[plt_index], nuc_count_name)
                        x_tickes = np.arange(0, 1.05, 0.05)
                        ax[plt_index].set_xticks(x_tickes)
                        
                        # result_list.append(ratio_attr) 
                    elif attr == archiveType.LMV:
                        new_attr:EV = getattr(design.investigator.investigator_results.lmv_values.lmv_comps[plt_index], nuc_count_name)
                        new_attr_value = new_attr.ev_normalized
                        # result_list.append(new_attr.ev_normalized)
                    
                    elif attr == archiveType.STATIC_PRIMES:
                        static_detector:StaticSystemDetector = StaticSystemDetector()
                        static_primes_nuc_count:PrimeNucCounts = static_detector.find_3prime_5prime_static_system(unbound_structure=design.investigator.lmv_references.weighted_structures.structs[plt_index],#design.investigator.lmv_references.mfe_structure,
                                                                                                                    bound_structure=source_data[index].fmn_folded_weighted) #design.investigator.lmv_references.weighted_structures.structs[plt_index])
                        static__nuc_ratio:float = float(getattr(static_primes_nuc_count, nuc_count_name)) / design.investigator.lmv_references.mfe_structure.nuc_count
                        new_attr_value = static__nuc_ratio
                        x_tickes = np.arange(0, 1.05, 0.05)
                        ax[plt_index].set_xticks(x_tickes)
                    
                    elif attr == archiveType.SNARE:
                        detector:MolecularSnareDetector = MolecularSnareDetector()
                        snare_result:SnareResults = detector.find_prime_moleculare_snare(moleculte_binding_sequence=snare_binding,
                                                             unbound_secondary_structure=design.investigator.lmv_references.weighted_structures.structs[plt_index],
                                                             bound_secondary_structure=source_data[index].fmn_folded_weighted)
                        if snare_result.snare_dectected is True:
                            snare_nuc_ratio:float = float(snare_result.snare_list[0].snare_stem_nuc_count) / design.investigator.lmv_references.mfe_structure.nuc_count
                            new_attr_value = snare_nuc_ratio
                            x_tickes = np.arange(0, 1.05, 0.05)
                            ax[plt_index].set_xticks(x_tickes)
                            
                  
                    
                    if design.investigator.number_structures[0] <= low_structs_num:
                        low_structs.append(new_attr_value)
                        result_fold_change_low_structs.append(design.design_info.wetlab_results.FoldChange)
                        baseline_low_structs.append(design.design_info.wetlab_results.Baseline_Subscore)
                        switch_low_structs.append(design.design_info.wetlab_results.Switch_Subscore)
                        folding_low_structs.append(design.design_info.wetlab_results.Folding_Subscore)
                        kdoff_low_structs.append(design.design_info.wetlab_results.KDOFF)
                        kdone_low_structs.append(design.design_info.wetlab_results.KDON)
                        eterna_low_structs.append(design.design_info.wetlab_results.Eterna_Score)
                        serena_basic_low_structs.append(design.investigator.basic_scores[0].total_score)
                        serena_advanced_low_structs.append(design.investigator.basic_scores[0].total_score + design.investigator.advanced_scores[0].total_score)
                        
                    elif design.investigator.number_structures[0] > low_structs_num and design.investigator.number_structures[0] <= med_structs_num:
                        med_structs.append(new_attr_value)
                        result_fold_change_med_structs.append(design.design_info.wetlab_results.FoldChange)
                        baseline_med_structs.append(design.design_info.wetlab_results.Baseline_Subscore)
                        switch_med_structs.append(design.design_info.wetlab_results.Switch_Subscore)
                        folding_med_structs.append(design.design_info.wetlab_results.Folding_Subscore)
                        kdoff_med_structs.append(design.design_info.wetlab_results.KDOFF)
                        kdone_med_structs.append(design.design_info.wetlab_results.KDON)
                        eterna_med_structs.append(design.design_info.wetlab_results.Eterna_Score)
                        serena_basic_med_structs.append(design.investigator.basic_scores[0].total_score)
                        serena_advanced_med_structs.append(design.investigator.basic_scores[0].total_score + design.investigator.advanced_scores[0].total_score)
                        
                    elif design.investigator.number_structures[0] > med_structs_num and design.investigator.number_structures[0] <= high_structs_num:
                        high_structs.append(new_attr_value)
                        result_fold_change_high_structs.append(design.design_info.wetlab_results.FoldChange)
                        baseline_high_structs.append(design.design_info.wetlab_results.Baseline_Subscore)
                        switch_high_structs.append(design.design_info.wetlab_results.Switch_Subscore)
                        folding_high_structs.append(design.design_info.wetlab_results.Folding_Subscore)
                        kdoff_high_structs.append(design.design_info.wetlab_results.KDOFF)
                        kdone_high_structs.append(design.design_info.wetlab_results.KDON)
                        eterna_high_structs.append(design.design_info.wetlab_results.Eterna_Score)
                        serena_basic_high_structs.append(design.investigator.basic_scores[0].total_score)
                        serena_advanced_high_structs.append(design.investigator.basic_scores[0].total_score + design.investigator.advanced_scores[0].total_score)
                        
                    elif design.investigator.number_structures[0] > high_structs_num:
                        obsurd_structs.append(new_attr_value)
                        result_fold_change_obsurd_structs.append(design.design_info.wetlab_results.FoldChange)
                        baseline_obsurd_structs.append(design.design_info.wetlab_results.Baseline_Subscore)
                        switch_obsurd_structs.append(design.design_info.wetlab_results.Switch_Subscore)
                        folding_obsurd_structs.append(design.design_info.wetlab_results.Folding_Subscore)
                        kdoff_obsurd_structs.append(design.design_info.wetlab_results.KDOFF)
                        kdone_obsurd_structs.append(design.design_info.wetlab_results.KDON)
                        eterna_obsurd_structs.append(design.design_info.wetlab_results.Eterna_Score)
                        serena_basic_obsurd_structs.append(design.investigator.basic_scores[0].total_score)
                        serena_advanced_obsurd_structs.append(design.investigator.basic_scores[0].total_score + design.investigator.advanced_scores[0].total_score)
                    # elif attr == archiveType.LMV_MFE:
                    #     good_list.append(getattr(design.investigator.investigator_results.lmv_values.lmv_comps[plt_index], nuc_count_name).ev_normalized)
                    # elif attr == archiveType.LMV_COMP:
                    #     good_list.append(getattr(design.investigator.investigator_results.lmv_values.lmv_comps[plt_index], nuc_count_name).ev_normalized)    
                    # result_fold_change.append(design.design_info.wetlab_results.FoldChange)
                    # baseline_scores.append(design.design_info.wetlab_results.Baseline_Subscore)
                    # switch_scores.append(design.design_info.wetlab_results.Switch_Subscore)
                    # folding_scores.append(design.design_info.wetlab_results.Folding_Subscore)
                    # kdoff_values.append(design.design_info.wetlab_results.KDOFF)
                    # kdone_values.append(design.design_info.wetlab_results.KDON)
                    # eterna_values.append(design.design_info.wetlab_results.Eterna_Score)
                    # serena_basic_values.append(design.investigator.basic_scores[0].total_score)
                    # serena_advanced_values.append(design.investigator.basic_scores[0].total_score + design.investigator.advanced_scores[0].total_score)
            
            filename_type:str = ''   
            if training is True:
                ax[plt_index].scatter(good_list, good_fold_change, color='green')
                ax[plt_index].scatter(bad_list, bad_change, c='red')
            else: 
                
                if score_type == ScoreType.FOLDCHANGE:               
                    # ax[plt_index].scatter(result_list, result_fold_change, c='blue')
                    ax[plt_index].scatter(low_structs, result_fold_change_low_structs, color='green', marker=MarkerStyle('o', 'none'), label=f'less than {str(low_structs_num)} stucts')
                    ax[plt_index].scatter(med_structs, result_fold_change_med_structs, color='blue',marker=MarkerStyle('o', 'none'),label=f'between {str(low_structs_num)} and {str(med_structs_num)} stucts')
                    ax[plt_index].scatter(high_structs, result_fold_change_high_structs, color='black' ,marker=MarkerStyle('o', 'none'),label=f'between {str(med_structs_num)} and {str(high_structs_num)} stucts')
                    ax[plt_index].scatter(obsurd_structs, result_fold_change_obsurd_structs, color='red' ,marker=MarkerStyle('o', 'none'),label=f'greater than {str(high_structs_num)} stucts')
                    ax[plt_index].set_ylabel("Foldchange")
                    filename_type = 'Foldchange'
                elif score_type == ScoreType.BASELINE:
                    # ax[plt_index].scatter(result_list, baseline_scores, color='red')
                    ax[plt_index].scatter(low_structs, baseline_low_structs, color='green',marker=MarkerStyle('o', 'none'), label=f'less than {str(low_structs_num)} stucts')
                    ax[plt_index].scatter(med_structs, baseline_med_structs, color='blue',marker=MarkerStyle('o', 'none'),label=f'between {str(low_structs_num)} and {str(med_structs_num)} stucts')
                    ax[plt_index].scatter(high_structs, baseline_high_structs, color='black',marker=MarkerStyle('o', 'none'),label=f'between {str(med_structs_num)} and {str(high_structs_num)} stucts')
                    ax[plt_index].scatter(obsurd_structs, baseline_obsurd_structs, color='red',marker=MarkerStyle('o', 'none'),label=f'greater than {str(high_structs_num)} stucts')
                    ax[plt_index].set_ylabel("Eterna Baseline subscore")
                    filename_type = 'Baseline'
                elif score_type == ScoreType.FOLDING:
                    # ax[plt_index].scatter(result_list, folding_scores, color='green')
                    ax[plt_index].scatter(low_structs, folding_low_structs, color='green',marker=MarkerStyle('o', 'none'), label=f'less than {str(low_structs_num)} stucts')
                    ax[plt_index].scatter(med_structs, folding_med_structs, color='blue',marker=MarkerStyle('o', 'none'),label=f'between {str(low_structs_num)} and {str(med_structs_num)} stucts')
                    ax[plt_index].scatter(high_structs, folding_high_structs, color='black',marker=MarkerStyle('o', 'none'),label=f'between {str(med_structs_num)} and {str(high_structs_num)} stucts')
                    ax[plt_index].scatter(obsurd_structs, folding_obsurd_structs, color='red',marker=MarkerStyle('o', 'none'),label=f'greater than {str(high_structs_num)} stucts')
                    ax[plt_index].set_ylabel("Eterna Folding subscore")
                    filename_type = 'Folding'
                elif score_type == ScoreType.SWITCH:
                    # ax[plt_index].scatter(result_list, switch_scores, color='brown')
                    ax[plt_index].scatter(low_structs, switch_low_structs, color='green',marker=MarkerStyle('o', 'none'), label=f'less than {str(low_structs_num)} stucts')
                    ax[plt_index].scatter(med_structs, switch_med_structs, color='blue',marker=MarkerStyle('o', 'none'),label=f'between {str(low_structs_num)} and {str(med_structs_num)} stucts')
                    ax[plt_index].scatter(high_structs, switch_high_structs, color='black',marker=MarkerStyle('o', 'none'),label=f'between {str(med_structs_num)} and {str(high_structs_num)} stucts')
                    ax[plt_index].scatter(obsurd_structs, switch_obsurd_structs, color='red',marker=MarkerStyle('o', 'none'),label=f'greater than {str(high_structs_num)} stucts')
                    ax[plt_index].set_ylabel("Eterna Switching subscore")
                    filename_type = 'Switching'
                elif score_type == ScoreType.KDON:
                    # ax[plt_index].scatter(result_list, kdone_values, color='purple')
                    ax[plt_index].scatter(low_structs, kdone_low_structs, color='green',marker=MarkerStyle('o', 'none'), label=f'less than {str(low_structs_num)} stucts')
                    ax[plt_index].scatter(med_structs, kdone_med_structs, color='blue',marker=MarkerStyle('o', 'none'),label=f'between {str(low_structs_num)} and {str(med_structs_num)} stucts')
                    ax[plt_index].scatter(high_structs, kdone_high_structs, color='black',marker=MarkerStyle('o', 'none'),label=f'between {str(med_structs_num)} and {str(high_structs_num)} stucts')
                    ax[plt_index].scatter(obsurd_structs, kdone_obsurd_structs, color='red',marker=MarkerStyle('o', 'none'),label=f'greater than {str(high_structs_num)} stucts')
                    ax[plt_index].set_ylabel("KDON")
                    filename_type = 'KDON'
                elif score_type == ScoreType.KDOFF:
                    # ax[plt_index].scatter(result_list, kdoff_values, color='black')
                    ax[plt_index].scatter(low_structs, kdoff_low_structs, color='green',marker=MarkerStyle('o', 'none'), label=f'less than {str(low_structs_num)} stucts')
                    ax[plt_index].scatter(med_structs, kdone_med_structs, color='blue',marker=MarkerStyle('o', 'none'),label=f'between {str(low_structs_num)} and {str(med_structs_num)} stucts')
                    ax[plt_index].scatter(high_structs, kdoff_high_structs, color='black',marker=MarkerStyle('o', 'none'),label=f'between {str(med_structs_num)} and {str(high_structs_num)} stucts')
                    ax[plt_index].scatter(obsurd_structs, kdoff_obsurd_structs, color='red',marker=MarkerStyle('o', 'none'),label=f'greater than {str(high_structs_num)} stucts')
                    ax[plt_index].set_ylabel("KDOFF")
                    filename_type = 'KDOFF'
                elif score_type == ScoreType.ETERNA:
                    # ax[plt_index].scatter(result_list, eterna_values, color='black')
                    ax[plt_index].scatter(low_structs, eterna_low_structs, color='green',marker=MarkerStyle('o', 'none'), label=f'less than {str(low_structs_num)} stucts')
                    ax[plt_index].scatter(med_structs, eterna_med_structs, color='blue',marker=MarkerStyle('o', 'none'),label=f'between {str(low_structs_num)} and {str(med_structs_num)} stucts')
                    ax[plt_index].scatter(high_structs, eterna_high_structs, color='black',marker=MarkerStyle('o', 'none'), label=f'between {str(med_structs_num)} and {str(high_structs_num)} stucts')
                    ax[plt_index].scatter(obsurd_structs, eterna_obsurd_structs, color='red',marker=MarkerStyle('o', 'none'),label=f'greater than {str(high_structs_num)} stucts')
                    ax[plt_index].set_ylabel("Eterna Score")
                    filename_type = 'Eterna'
                elif score_type == ScoreType.BASIC:
                    # ax[plt_index].scatter(result_list, serena_basic_values, color='black')
                    ax[plt_index].scatter(low_structs, serena_basic_low_structs, color='green',marker=MarkerStyle('o', 'none'), label=f'less than {str(low_structs_num)} stucts')
                    ax[plt_index].scatter(med_structs, serena_basic_med_structs, color='blue',marker=MarkerStyle('o', 'none'),label=f'between {str(low_structs_num)} and {str(med_structs_num)} stucts')
                    ax[plt_index].scatter(high_structs, serena_basic_high_structs, color='black',marker=MarkerStyle('o', 'none'),label=f'between {str(med_structs_num)} and {str(high_structs_num)} stucts')
                    ax[plt_index].scatter(obsurd_structs, serena_basic_obsurd_structs, color='red',marker=MarkerStyle('o', 'none'),label=f'greater than {str(high_structs_num)} stucts')
                    ax[plt_index].set_ylabel("Serena Basic Score")
                    ax[plt_index].set_ylim(-20, 20)
                    filename_type = 'Basic'
                elif score_type == ScoreType.ADVANCED:
                    # ax[plt_index].scatter(result_list, serena_advanced_values, color='black')
                    ax[plt_index].scatter(low_structs, serena_advanced_low_structs, color='green',marker=MarkerStyle('o', 'none'), label=f'less than {str(low_structs_num)} stucts')
                    ax[plt_index].scatter(med_structs, serena_advanced_med_structs, color='blue',marker=MarkerStyle('o', 'none'),label=f'between {str(low_structs_num)} and {str(med_structs_num)} stucts')
                    ax[plt_index].scatter(high_structs, serena_advanced_high_structs, color='black',marker=MarkerStyle('o', 'none'),label=f'between {str(med_structs_num)} and {str(high_structs_num)} stucts')
                    ax[plt_index].scatter(obsurd_structs, serena_advanced_obsurd_structs, color='red',marker=MarkerStyle('o', 'none'),label=f'greater than {str(high_structs_num)} stucts')
                    ax[plt_index].set_ylabel("Serena Basic plus Advanced Score")
                    filename_type = 'Advanced' 
                    ax[plt_index].set_ylim(-20, 20)   
                
            if attr == archiveType.RATIO or attr == archiveType.STATIC_PRIMES:
                ax[plt_index].legend(loc="upper left")
            else:
                ax[plt_index].legend(loc="best")
            
            kcal_delta = plt_index +1
            stuff = f'{kcal_delta-1}kcal to {kcal_delta}kcal delta from MFE'
            ax[plt_index].set_xlabel(stuff)
            
            ax[plt_index].set_xlim(-.3, x_range+.05)
            
            
            
        # fig.tight_layout()
        plt.subplots_adjust(left=0.1,
                bottom=.1,  
                top=.9, 
                )       
        
        save_dir:str = f'/home/rnauser/repo/Serena-Local-Minima-Variation-Tool/src/tests/bin/{timestr}'
        if os.path.isdir(save_dir) == False:
            os.makedirs(save_dir)
        plt.savefig(f'{save_dir}/{nuc_count_name}_{filename_type}_{timestr}.png')
        plt.close()
        # plt.show()
        
        
def plot_investigator(sublab:str, test_name:str, cluster_size_threshold:int, pnas_path:Path, round:str, archive_path:Path, source_archive_path:Path, snare_binding:str = None ):
    # sublab:str = 'SSNG1'
    test_name:str = f'{test_name}_cluster_{cluster_size_threshold}'
    timestr:str = f'{sublab}_{test_name}_{time.strftime("%Y%m%d-%H%M%S")}'
    plot_investigaot:InvestigatorReportGeneration = InvestigatorReportGeneration()
    pnas:ProcessPNAS = ProcessPNAS()
    
    # archive_path:str = f'/home/rnauser/test_data/serena/R101_PNAS/computational_data/{sublab}' #'/home/rnauser/test_data/serena/computatation/computational_data'#
    # source_archive_path:str =  f'/home/rnauser/test_data/serena/R101_PNAS/raw_fold_data/rna95_nupack3/{sublab}'
    archive_path = archive_path.joinpath(sublab)
    source_archive_path = source_archive_path.joinpath(sublab)
    
    
    new_sara:Sara2API = Sara2API()
    puzzle_data: puzzleData
    pandas_sheet: DataFrame
    # pnas_path:str = '/home/rnauser/test_data/serena/R101_PNAS/source/pnas.2112979119.sd01.xlsx' #'/home/rnauser/test_data/pnas_testing_tweak.xlsx'#
    puzzle_data, pandas_sheet = new_sara.ProcessLab(path=pnas_path.as_posix(),
                                                    designRound_sheet=round,#"Round 7 (R101)", #'R101 Filtered good bad',#'R101 Filtered good bad',#
                                                    sublab_name=''
                                                    )
    pnas_data:List[ArchiveInvestigatorData] = []
    source_data:List[ArchiveData] = []
    
    # ignor_list = [6387992, 6374121, 6374127, 6388849, 6388889, 6391463]
    dir_list = os.listdir(archive_path)
    flag = 0
    for design in puzzle_data.designsList: 
        
        if str(design.design_info.DesignID) in dir_list and design.wetlab_results.NumberOfClusters1 > cluster_size_threshold: 
        
            archived_data:ArchiveInvestigatorData = ArchiveInvestigatorData(design_info=design)
            
            archived_data = pnas.archive_investigation_computations(dest_folder=archive_path.as_posix(),
                                                                    flow=ArchiveFlow.GET,
                                                                    data=archived_data)
            
            temp_archive:ArchiveData = ArchiveData(design_info=design)
            # found_source_data:ArchiveData = pnas.archive_ensemble_data(dest_folder=source_archive_path,
            #                             flow=ArchiveFlow.GET,
            #                             data=temp_archive)
            
            backup_records:ArchiveSecondaryStructureList = ArchiveSecondaryStructureList(working_folder=source_archive_path.as_posix(),
                                             var_name=str(design.design_info.DesignID),
                                             use_db=True)
            temp_archive.fmn_folded_weighted = backup_records.data.fmn_folded_weighted
            
            if archived_data.design_info.wetlab_results.Eterna_Score == 100:
            
                source_data.append(temp_archive)
                pnas_data.append(archived_data)
                break
        
            # if flag > 5:
            #     break
            # flag +=1
    
    ratio_value:int = 1
    
    # if archived_data.design_info.wetlab_results.NumberOfClusters1 > 0:#200:
    snare_test_name:str = "moleculare_snare_nuc_ratio"
    if snare_binding != None:
        for enumerator in ScoreType:
            plot_investigaot.generate_nuc_count_plot(x_range=ratio_value,
                                                    data=pnas_data,
                                                    source_data=source_data,
                                                    attr=archiveType.SNARE, 
                                                    nuc_count_name=snare_test_name, x_string="Ratio of molecular snare static nucs to total nucs",
                                                    training=False,
                                                    score_type=enumerator,
                                                    timestr=timestr,
                                                    snare_binding=snare_binding)
    
    for item in ['static_stem_nuc_count', 'static_loop_nuc_count', 'prime_static_nuc_count_total']:
        for enumerator in ScoreType:
            plot_investigaot.generate_nuc_count_plot(x_range=ratio_value,
                                                    data=pnas_data,
                                                    source_data=source_data,
                                                    attr=archiveType.STATIC_PRIMES, 
                                                    nuc_count_name=item, x_string="Ratio of static 5' and 3' static nucleotides",
                                                    training=False,
                                                    score_type=enumerator,
                                                    timestr=timestr)

    for count in  archived_data.investigator.investigator_results.comp_nuc_counts.comparison_nuc_counts[0].__dict__:
        
        for enumerator in ScoreType:
        
            plot_investigaot.generate_nuc_count_plot(x_range=archived_data.investigator.investigator_results.comp_nuc_counts.comparison_nuc_counts[0].num_nucs,
                                                data=pnas_data,
                                                source_data=source_data, 
                                                attr=archiveType.COUNT, 
                                                nuc_count_name=count, x_string="Count of Nucleotides",
                                                training=False,
                                                score_type=enumerator,
                                                timestr=timestr)
     
    
    for ratio in  archived_data.investigator.investigator_results.comparison_eval_results.ratios[0].__dict__:
        
        for enumerator in ScoreType:
        
            plot_investigaot.generate_nuc_count_plot(x_range=ratio_value,
                                                data=pnas_data,
                                                source_data=source_data, 
                                                attr=archiveType.RATIO, 
                                                nuc_count_name=ratio, x_string="Ratio of nucleotide position counts",
                                                training=False,
                                                score_type=enumerator,
                                                timestr=timestr)
      
    for lmv_rel in archived_data.investigator.investigator_results.lmv_values.lmv_comps[0].__dict__:
        
        for enumerator in ScoreType:
            plot_investigaot.generate_nuc_count_plot(x_range=40,
                                                data=pnas_data,
                                                source_data=source_data, 
                                                attr=archiveType.LMV, 
                                                nuc_count_name=lmv_rel, x_string="LMV of group",
                                                training=False,
                                                score_type=enumerator,
                                                timestr=timestr)
                
  

def run_plot_investigator():
    parser = argparse.ArgumentParser(description='Get and process R101 PNAS data')
    
    parser.add_argument('--pnas', 
                        type=Path,
                        required=True,
                        help='Path to the pnas file for analysis')
    
    parser.add_argument('--round', 
                        type=str,
                        default='Round 7 (R101)',
                        required=False,
                        help='Round to run')
    
    parser.add_argument('--sublab', 
                        type=str,
                        default='',
                        required=True,
                        help='sublab in R101 to generate and archive enemble data for')
    
    parser.add_argument('--test-name', 
                        type=str,
                        default='',
                        required=True,
                        help='name for the test')
    
    parser.add_argument('--source-data', 
                        type=Path,
                        required=True,
                        help='Path to the data nut squirrel archive folder')
    
    parser.add_argument('--computations', 
                        type=Path,
                        required=True,
                        help='Path to computationse folder')
    
    parser.add_argument('--do-weighted', 
                        action="store_true",
                        required=False,
                        help='Use weighted struct from Vienna2 enemble')
    
    parser.add_argument('--cluster', 
                        type=int,
                        default=500000,
                        required=False,
                        help='sublab in R101 to generate and archive enemble data for')
    
    parser.add_argument('--snare-binding', 
                        type=str,
                        default=None,
                        required=False,
                        help='moleculare snare binding sequence')
    
    
    
    args = parser.parse_args()
    

    
    
    # print(args.do_agressive)
    
    plot_investigator(sublab=args.sublab,
                      test_name=args.test_name,
                      cluster_size_threshold=args.cluster,
                      pnas_path=args.pnas,
                      round=args.round,
                      archive_path=args.computations,
                      source_archive_path=args.source_data,
                      snare_binding=args.snare_binding)
           
        
plot_investigator(sublab='SSNG1',
                  test_name='molecular_snare',
                  cluster_size_threshold=200,
                  pnas_path=Path('/home/rnauser/test_data/serena/R101_PNAS/source/pnas.2112979119.sd01.xlsx'),
                  round='Round 7 (R101)',
                  archive_path=Path('/home/rnauser/test_data/serena/R101_PNAS/computational_data/'),
                  source_archive_path=Path('/home/rnauser/test_data/serena/R101_PNAS/raw_fold_data/rna95_nupack3/'),
                  snare_binding='AGGAUAUAGAAGG')