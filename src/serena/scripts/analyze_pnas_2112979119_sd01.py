"""
File for managing analysis of dataset pnas.2112979119.sd01 from https://www.pnas.org/doi/full/10.1073/pnas.2112979119
"""

from enum import Enum
from typing import List
from pathlib import Path
from serena.utilities.ensemble_structures import Sara2SecondaryStructure, Sara2StructureList

from serena.bin.backup_serena import ArchiveSecondaryStructureList

class ArchiveFlow(Enum):
    GET="GET"
    PUT="PUT"


class ProcessPNAS():
    
    def __init__(self) -> None:
        pass
    
    def get_nupack_ensemble_structs(self, pnas_dataset_path:Path, round:str, sublab:str)->List[Sara2StructureList]:
        pass
    
    def archive_ensemble_structs(self, structs:List[Sara2StructureList], dest_folder:Path, flow:ArchiveFlow):
        pass
    
    def switchyness_analysis(self, designs_structures:List[Sara2StructureList]):
        pass