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

from enum import Enum

import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import matplotlib

from serena.utilities.ensemble_variation import EnsembleVariation, EVResult

from serena.analysis.investigator import (
    InvestigatorResults,
    RatioResults,
    
)

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
    ComparisonLMVResponse
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

class archiveType(Enum):
    RATIO='RATIO'
    COUNT="COUNT"

class InvestigatorReportGeneration():
    
    def __init__(self) -> None:
        pass
    
            
    def generate_nuc_count_plot(self, data:List[ArchiveInvestigatorData], attr:archiveType, nuc_count_name:str, x_string:str, x_range:float):
        timestr = time.strftime("%Y%m%d-%H%M%S")
        
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
            
            for design in data:
                if design.design_info.design_info.Puzzle_Name == 'good':
                    if attr == archiveType.COUNT:
                        good_list.append(getattr(design.investigator.investigator_results.comp_nuc_counts.comparison_nuc_counts[plt_index], nuc_count_name)) 
                    elif attr == archiveType.RATIO:
                        good_list.append(getattr(design.investigator.investigator_results.comparison_eval_results.ratios[plt_index], nuc_count_name)) 
                    good_fold_change.append(design.design_info.wetlab_results.FoldChange)
                elif design.design_info.design_info.Puzzle_Name == 'bad':
                    if attr == archiveType.COUNT:
                        bad_list.append(getattr(design.investigator.investigator_results.comp_nuc_counts.comparison_nuc_counts[plt_index], nuc_count_name)) 
                    elif attr == archiveType.RATIO:
                        bad_list.append(getattr(design.investigator.investigator_results.comparison_eval_results.ratios[plt_index], nuc_count_name))  
                    bad_change.append(design.design_info.wetlab_results.FoldChange)
            
            ax[plt_index].scatter(good_list, good_fold_change, c='green')
            ax[plt_index].scatter(bad_list, bad_change, c='red')
            kcal_delta = plt_index +1
            stuff = f'{kcal_delta}kcal delta from MFE'
            ax[plt_index].set_xlabel(stuff)
            ax[plt_index].set_ylabel("Fold Change")
            ax[plt_index].set_xlim(0, x_range)
            
        # fig.tight_layout()
        plt.subplots_adjust(left=0.1,
                bottom=.1,  
                top=.9, 
                )       
        
        plt.savefig(f'/home/rnauser/repo/Serena-Local-Minima-Variation-Tool/src/tests/bin/{nuc_count_name}_{timestr}.png')
        # plt.show()
        
        
def plot_investigator():
    
    plot_investigaot:InvestigatorReportGeneration = InvestigatorReportGeneration()
    pnas:ProcessPNAS = ProcessPNAS()
    
    archive_path:str = '/home/rnauser/test_data/serena/computatation/computational_data'
    
    new_sara:Sara2API = Sara2API()
    puzzle_data: puzzleData
    pandas_sheet: DataFrame
    pnas_path:str = '/home/rnauser/test_data/pnas_testing_tweak.xlsx'
    puzzle_data, pandas_sheet = new_sara.ProcessLab(path=pnas_path,
                                                    designRound_sheet="R101 Filtered good bad",
                                                    sublab_name=''
                                                    )
    pnas_data:List[ArchiveInvestigatorData] = []
    
    ignor_list = [6387992, 6374121, 6374127, 6388849, 6388889, 6391463]
    dir_list = os.listdir(archive_path)
    flag = 0
    for design in puzzle_data.designsList: 
        
        if str(design.design_info.DesignID) in dir_list and design.design_info.DesignID not in ignor_list:   
        
            archived_data:ArchiveInvestigatorData = ArchiveInvestigatorData(design_info=design)
            
            archived_data = pnas.archive_investigation_computations(dest_folder=archive_path,
                                                                    flow=ArchiveFlow.GET,
                                                                    data=archived_data)
            pnas_data.append(archived_data)
            # time.sleep(1)
            # flag += 1
            # if flag > 5:
            #     break
    
    for count in  archived_data.investigator.investigator_results.comp_nuc_counts.comparison_nuc_counts[0].__dict__:
        
        plot_investigaot.generate_nuc_count_plot(x_range=archived_data.investigator.investigator_results.comp_nuc_counts.comparison_nuc_counts[0].num_nucs,
                                                 data=pnas_data, 
                                                 attr=archiveType.COUNT, 
                                                 nuc_count_name=count, x_string="Count of Nucleotides")
    
    for ratio in  archived_data.investigator.investigator_results.comparison_eval_results.ratios[0].__dict__:
        
        plot_investigaot.generate_nuc_count_plot(x_range=1,
                                                 data=pnas_data, 
                                                 attr=archiveType.RATIO, 
                                                 nuc_count_name=ratio, x_string="Ratio of nucleotide position counts")
        
plot_investigator()