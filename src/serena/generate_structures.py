"""
This file is the entry point for constructing serena/sara structures
"""

from typing import List

from serena.utilities.ensemble_structures import Sara2StructureList, Sara2SecondaryStructure
from serena.utilities.ensemble_groups import SingleEnsembleGroup, MultipleEnsembleGroups, EnsembleSwitchStateMFEStructs

class SecondaryStructures():
    """
    Class to genereate the secondary structure
    framework used by serena and sara
    """

    def make_secondary_structure(primary_structure:str, secondary_structure:str, free_energy:float, stack_free_energy:float)->Sara2SecondaryStructure:
        return Sara2SecondaryStructure(sequence=primary_structure,
                                       structure=secondary_structure,
                                       freeEnergy=free_energy,
                                       stackEnergy=stack_free_energy
                                       )
    
    def make_secondary_strucuture_list(secondary_structures_list: List[Sara2SecondaryStructure])->Sara2StructureList:
        structure_list:Sara2StructureList = Sara2StructureList()
        for structure in secondary_structures_list:
            structure_list.add_structure(structure)
        return Sara2StructureList
    
class EnsembleGroups():
    """
    Class for generating the ensemble groups consumed by serena and sara
    """

    def make_switch_mfe_states_from_secondary_strucures(switched_state_mfe_structure:Sara2SecondaryStructure, non_switch_mfe_structure:Sara2SecondaryStructure):
        return EnsembleSwitchStateMFEStructs(non_switch_mfe_struct=non_switch_mfe_structure,
                                             switched_mfe_struct=switched_state_mfe_structure)
    
    def make_singel_ensemble_group(ensemble_structures:Sara2StructureList, mfe_switch_structures:EnsembleSwitchStateMFEStructs, kcal_start:float, kcal_end:float):
        single_ensemble_group:SingleEnsembleGroup = SingleEnsembleGroup()
        single_ensemble_group.group = ensemble_structures
        single_ensemble_group.switch_state_structures = mfe_switch_structures
        single_ensemble_group.kcal_start = kcal_start
        single_ensemble_group.kcal_end = kcal_end
        single_ensemble_group.kcal_span = kcal_end - kcal_start
        return single_ensemble_group
    
    def make_multiple_ensemple_groups(ensemble_groups:List[SingleEnsembleGroup], mfe_switch_structures:EnsembleSwitchStateMFEStructs):
        multiple_ensemble_group:MultipleEnsembleGroups = MultipleEnsembleGroups(switch_state_structures=mfe_switch_structures)
        for group in SingleEnsembleGroup:
            multiple_ensemble_group.add_group(group=group)
        return multiple_ensemble_group