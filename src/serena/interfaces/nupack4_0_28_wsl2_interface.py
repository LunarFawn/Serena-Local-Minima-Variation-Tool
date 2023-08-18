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

from serena.utilities.ensemble_structures import  Sara2SecondaryStructure, Sara2StructureList, MakeSecondaryStructures
from serena.utilities.ensemble_groups import SingleEnsembleGroup, MultipleEnsembleGroups, MakeEnsembleGroups, EnsembleSwitchStateMFEStructs

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
    
    def load_nupack_subopt_as_ensemble(self, span_structures:Sara2StructureList, settings: NupackSettings, switch_state:EnsembleSwitchStateMFEStructs):
        make_ensemble: MakeEnsembleGroups = MakeEnsembleGroups()
        make_structs: MakeSecondaryStructures = MakeSecondaryStructures()
        mfe_energy:float =  span_structures.mfe_freeEnergy

        #this is for increments of 1 kcal need to do fraction
        num_groups: int = int(settings.kcal_span_from_mfe / settings.Kcal_unit_increments)
        remainder: int = settings.kcal_span_from_mfe % settings.Kcal_unit_increments

        groups_list : List[Sara2StructureList] = []
        groups_index_used: List[bool] = []
        groups_dict: Dict[int, Sara2StructureList] = {}
        group_values: List[float] = []

        #this fills up the list of energy deltas to publich EV's for
        current_energy: float = mfe_energy
        group_values.append(current_energy)
        for index in range(num_groups):
            current_energy = current_energy + settings.Kcal_unit_increments
            group_values.append(current_energy)
        
        #now initialize the groups_list
        for index in range(len(group_values)-1):
            group: Sara2StructureList = Sara2StructureList()
            groups_list.append(group)
            groups_index_used.append(False)
            groups_dict[index+1] = group

        num_sara_struct: int = span_structures.num_structures
        for sara_index in range(0,num_sara_struct):
            sara_structure: Sara2SecondaryStructure = span_structures.sara_stuctures[sara_index]
            current_energy = sara_structure.freeEnergy

            #need to do this because there are two indexes need to look at each 
            #loop and want to avoid triggering a list index overrun
            for group_index in range(len(group_values)-1):
                #remember we are dealing with neg kcal so its you want to 
                min_energy: float = group_values[group_index]
                max_energy: float = group_values[group_index+1]
                if current_energy >= min_energy and current_energy < max_energy:
                    groups_list[group_index].add_structure(sara_structure)
                    groups_index_used[group_index] = True            
    
        single_groups: List[SingleEnsembleGroup] = []
        
        for group_index in range(len(groups_list)):
            group = groups_list[group_index]
            start_value = group_values[group_index] - settings.Kcal_unit_increments
            end_value = group_values[group_index]
            this_group:SingleEnsembleGroup = make_ensemble.make_singel_ensemble_group(ensemble_structures=group,
                                                                                      mfe_switch_structures=switch_state,
                                                                                      kcal_start=start_value,
                                                                                      kcal_end=end_value)
            single_groups.append(this_group)

        ensemble_groups: MultipleEnsembleGroups = make_ensemble.make_multiple_ensemple_groups(ensemble_groups=single_groups,
                                                                                              mfe_switch_structures=switch_state)
        
        return ensemble_groups
            
        