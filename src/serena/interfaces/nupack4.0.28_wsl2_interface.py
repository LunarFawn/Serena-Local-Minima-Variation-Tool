"""
File that provides an interface to standard python nupack4 install on linux
Intention is to be able to access this via a docker that will be accessable on docker hub
that has nupack4 setup and ready for this project to consume
"""
from typing import List, Dict
from dataclasses import dataclass
from datetime import datetime

from enum import Enum
from nupack import *

from serena.utilities.ensemble_structures import  Sara2SecondaryStructure, Sara2StructureList
from serena.utilities.ensemble_groups import SingleEnsembleGroup, MultipleEnsembleGroups

class MaterialParameter(Enum):
    NONE = 0
    #"Based on [Mathews99] and [Lu06] with additional parameters [Xia98,Zuker03] including coaxial stacking [Mathews99,Turner10] and dangle stacking [Serra95,Zuker03,Turner10] in 1M Na+."
    rna06_nupack4 = 1
    #"Based on [Serra95] with additional parameters [Zuker03] including coaxial stacking [Mathews99,Turner10] and dangle stacking [Serra95,Zuker03,Turner10] in 1M Na+."
    rna95_nupack4 = 2    
    #"Parameters from [Mathews99] with terminal mismatch free energies in exterior loops and multiloops replaced by two dangle stacking free energies. Parameters are provided only for 37 ∘C."
    rna99_nupack3 = 3    
    #"Same as rna95 except that terminal mismatch free energies in exterior loops and multiloops are replaced by two dangle stacking free energies."
    rna95_nupack3 = 4

@dataclass
class NupackSettings():  
    material_param:MaterialParameter = MaterialParameter.NONE
    temp_C: int = 0
    kcal_span_from_mfe:int = 0
    Kcal_unit_increments: float = 0
    sequence:str = ''

class NUPACK4Interface():
    """
    Class for nupack4 interface for sara2 logic intended for serena package
    """

    def __init__(self) -> None:
        pass

    def select_material_parameters(self, parameters:MaterialParameter):
        param:str = ''        
        if parameters == MaterialParameter.rna06_nupack4:
            param = "rna06"
        elif parameters == MaterialParameter.rna95_nupack4:
            param = "rna95"
        elif parameters == MaterialParameter.rna99_nupack3:
            param = "rna99-nupack3"
        elif parameters == MaterialParameter.rna95_nupack3:
            param = "rna95-nupack3"            
        return param                                                      

    def select_model(self, material_param:MaterialParameter, temp_C:int):
        param: str = self.select_material_parameters(material_param)
        my_model = Model(material=param, celsius=temp_C)
        return my_model
    
    def get_subopt_energy_gap(self, material_param:MaterialParameter, temp_C:int, sequence_string:str, energy_delta_from_MFE: int):
        #run through subopt
        my_model = self.select_model(material_param, temp_C)
        kcal_group_structures_list: Sara2StructureList = Sara2StructureList()
        ensemble_kcal_group= subopt(strands=sequence_string, model=my_model, energy_gap=energy_delta_from_MFE)

        for i,kcal_group_elementInfo in enumerate(ensemble_kcal_group):
                  
            #get all the structures and energis pulled and prepped for proccessing and add them tot eh dict and the list               
            structure = str(kcal_group_elementInfo.structure)
            freeEnergy = float(kcal_group_elementInfo.energy)
            stackEnergy = float(kcal_group_elementInfo.stack_energy)
            structure_info: Sara2SecondaryStructure = Sara2SecondaryStructure(sequence=sequence_string, structure=structure, 
                                                                              freeEnergy=freeEnergy, stackEnergy=stackEnergy)
            kcal_group_structures_list.add_structure(structure_info)

        return kcal_group_structures_list