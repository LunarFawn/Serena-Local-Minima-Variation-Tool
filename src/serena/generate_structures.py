"""
This file is the entry point for constructing serena/sara structures
"""

from typing import List

from serena.utilities.ensemble_structures import Sara2StructureList, Sara2SecondaryStructure

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
            structure_list.add_structure(secondary_structures_list)
        return Sara2StructureList