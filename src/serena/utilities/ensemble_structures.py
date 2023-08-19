"""
Sara2 api for accessing and manipulating secondary structures 
in dot parenthisis form
copyright 2023 GrizzlyEngineer
"""
from typing import List, Dict, NamedTuple
import struct
import pandas as pd
import sys
import openpyxl
from copy import deepcopy
from dataclasses import dataclass
from datetime import datetime, timedelta
import threading
import time
from collections import namedtuple

@dataclass
class KcalRanges():
    """
    Class to hold kcal ranges
    """
    start: float = 0
    stop: float = 0

class Sara2SecondaryStructure(object):
    """
    Sara 2 Secondary Structure that is used to hold all the info for
    each secondary structure in ensemble
    """

    def __init__(self, sequence:str = '', structure: str = '', freeEnergy: float = 0, stackEnergy: float = 0) -> None:
        self._sequence: str = sequence
        self._structure: str = structure
        self._freeEnergy: float = freeEnergy
        self._stackEnergy: float = stackEnergy
        #self._nuc_count: int = len(sequence)

    @property
    def sequence(self):
        """
        Returns the sequence as a string
        """
        return self._sequence

    @sequence.setter
    def sequence(self, primary_struc: str):
        """
        Sets the sequence using string
        """
        self._sequence = primary_struc

    @property
    def structure(self):
        """
        Returns the secondary strucuture in dot parens notation
        """
        return self._structure

    @structure.setter
    def structure(self, dot_parens: str):
        """
        Sets the secondary structure using dot parense notation string
        """
        self._structure = dot_parens

    @property
    def freeEnergy(self):
        """
        Returns the total free energy as float
        """
        return self._freeEnergy

    @freeEnergy.setter
    def freeEnergy(self, energy: float):
        """
        Sets the total free energy with float
        """
        self._freeEnergy = energy

    @property
    def stackEnergy(self):
        """
        Returns the stack energy as float
        """
        return self._stackEnergy

    @stackEnergy.setter
    def stackEnergy(self, energy: float):
        """
        Sets the stack energy with float
        """
        self._stackEnergy = energy

    @property
    def nuc_count(self):
        """
        Returns the number of nucleotides as a int
        """
        return len(self._sequence)


class Sara2StructureList(object):
    """
    Sara2 Structure List that holds all the Sar2SecondaryStructurs
    that represent the ensemble in raw form
    """
    def __init__(self) -> None:
        self._sara_structures_list: List[Sara2SecondaryStructure] = []
        self._structures: List[str] = []
        self._freeEnergy_list: list[float] = []
        self._stackEnergy_list: list[float] = []
        self._min_freeEnergy: float = 0
        self._max_freeEnergy: float = 0
        self._min_stackEnergy: float = 0
        self._max_stackEnergy: float = 0
        self._num_structures: int = 0
        #self._nuc_count: int = 0
        #self._mfe_structure: str = ''
        #self._mfe_freeEnergy: float = 0
        #self._mfe_stackEnergy: float = 0
        self._freeEnergy_span:float = 0
        self._stackEnergy_span:float = 0
        self._weighted_structure:str = ''

    def process_energy(self):
        """
        Process min and max energies in list as well
        as populate counts. It always ran after adding 
        structure
        """
            #now populate min and max
        #do free energy
        if len(self._freeEnergy_list) == 0:
            self._min_freeEnergy = 0
            self._max_freeEnergy = 0
        else:
            self._min_freeEnergy = min(self._freeEnergy_list)
            self._max_freeEnergy = max(self._freeEnergy_list)

        self._freeEnergy_span = self._max_freeEnergy - self._min_freeEnergy
        #do stack energy

        if len(self._stackEnergy_list) == 0:
            self._min_stackEnergy = 0
            self._max_stackEnergy = 0
        else:
            self._min_stackEnergy = min(self._stackEnergy_list)
            self._max_stackEnergy = max(self._stackEnergy_list)
        self._stackEnergy_span = self._max_stackEnergy - self._min_stackEnergy

        #now count
        self._num_structures = len(self._sara_structures_list)

    def add_structure(self, structure: Sara2SecondaryStructure):
        """
        main way to add a structure to the list
        """
        self._sara_structures_list.append(structure)
        #self._structures.append(structure.structure)
        self._freeEnergy_list.append(structure.freeEnergy)
        self._stackEnergy_list.append(structure.stackEnergy)
        self.process_energy()

    def remove_structure(self, index:int):
        """
        remove a structure from memory
        """
        del self._structures[index]
        del self._freeEnergy_list[index]
        del self._stackEnergy_list[index]
        self.process_energy()            

    @property
    def mfe_structure(self):
        """
        Returns the mfe secibdary structure as a string
        """
        structure:str = ''
        if len(self.sara_stuctures) > 0:
           structure = self.sara_stuctures[0].structure
        return structure 

    @property
    def mfe_freeEnergy(self):
        """
        Returns the mfe total free energy as float
        """
        energy: float = 0
        if len(self.sara_stuctures) > 0:
            energy = self.sara_stuctures[0].freeEnergy
        return energy

    @property
    def mfe_stackEnergy(self):
        """
        Returns the mfe stack energy as float
        """
        energy: float = 0
        if len(self.sara_stuctures) > 0:
            energy = self.sara_stuctures[0].stackEnergy
        return energy

    @property
    def nuc_count(self):
        """
        Returns the total number of nucleotides as int
        """
        count: int = 0
        if len(self.sara_stuctures) > 0:
            count = self.sara_stuctures[0].nuc_count  
        return count

    @property
    def sara_stuctures(self):
        """
        Returns the sara structures that make up list
        """
        return self._sara_structures_list

    @sara_stuctures.setter   
    def sara_stuctures(self, structs_list: List[Sara2SecondaryStructure]):
        """
        Sets the sara structures list using a List of Sara2Structures
        """
        #reset list
        self._sara_structures_list=[]
        #fill it in now
        for struc in structs_list:
            self.add_structure(struc)

    @property
    def max_free_energy(self):
        """
        Returns the maximum free energy of the structures in the list
        """
        return self._max_freeEnergy

    @property
    def min_free_energy(self):
        """
        Returns the minimum free energy of the structures in the list
        """
        return self._min_freeEnergy

    @property
    def max_stack_energy(self):
        """
        Returns the maximum stack energy of the structures in the list
        """
        return self._max_stackEnergy

    @property
    def min_stack_energy(self):
        """
        Returns the minimum stack energy of the structures in the list
        """
        return self._min_stackEnergy

    @property
    def num_structures(self):
        """
        Returns the number of structures in the list
        """
        return self._num_structures

    @property
    def freeEnergy_span(self):
        """
        Returns the span of the free energy of the structures in the list
        """
        return self._freeEnergy_span

    @property
    def stackEnergy_span(self):
        """
        Returns the span of the stack energy of the structures in the list
        """
        return self._stackEnergy_span 

    @property
    def weighted_structure(self):
        """
        Returns the weighted structure as a string
        """
        return self._weighted_structure

    @weighted_structure.setter        
    def weighted_structure(self, structure: str):
        """
        sets the weigthed structure
        """
        self._weighted_structure = structure

class MakeSecondaryStructures():
    """
    Class to genereate the secondary structure
    framework used by serena and sara
    """
    def make_secondary_structure(self, primary_structure:str, secondary_structure:str, free_energy:float, stack_free_energy:float)->Sara2SecondaryStructure:
        """
        Function to make a secondary structue
        """
        return Sara2SecondaryStructure(sequence=primary_structure,
                                       structure=secondary_structure,
                                       freeEnergy=free_energy,
                                       stackEnergy=stack_free_energy
                                       )

    def make_secondary_strucuture_list(self, secondary_structures_list: List[Sara2SecondaryStructure])->Sara2StructureList:
        """
        Function to make a secondary structure list
        """
        structure_list:Sara2StructureList = Sara2StructureList()
        for structure in secondary_structures_list:
            structure_list.add_structure(structure)
        return structure_list
