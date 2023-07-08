
from dataclasses import dataclass
from typing import List
import collections

from serena.utilities.ensemble_groups import SingleEnsembleGroup
from serena.utilities.ensemble_structures import Sara2StructureList, Sara2SecondaryStructure


class WeightedStructure():

    def __init__(self,source_structures: Sara2StructureList ) -> None:        
        self._source_structures: Sara2StructureList = source_structures
        self._weighted_structure: Sara2SecondaryStructure = self.make_weighted_struct(source_structures)

    @property
    def weighted_structure(self):
        return self._weighted_structure
    
    @property
    def source_structures(self):
        return self._source_structures
    
    @source_structures.setter
    def source_structures(self, structures: Sara2StructureList):
        self._source_structures = structures
        self._weighted_structure: Sara2SecondaryStructure = self.make_weighted_struct(structures)

    def make_weighted_struct(self, structure_list: Sara2StructureList)->Sara2SecondaryStructure:
        is_bond_value: int = 2
        not_bond_value: int = -1

        nuc_poistion_values: List[int] = []
        nuc_pairs_comp_list: List[List[str]] = []
        good_nucs_each_pos: List[bool] = []

        struct_count: int = structure_list.num_structures

        for nucIndex in range(structure_list.nuc_count):
            nuc_poistion_values.append(0)
            pairs_list: List[str] = []            
            nuc_pairs_comp_list.append(pairs_list)
            #good_nucs_each_pos.append(False)

        for struct in structure_list.sara_stuctures:
            for nucIndex in range(structure_list.nuc_count):
                nuc_bond_type:str = struct.structure[nucIndex]
                nuc_pairs_comp_list[nucIndex].append(nuc_bond_type)
                adder: int = 0
                if nuc_bond_type == '.':
                    adder = not_bond_value
                else:
                    adder = is_bond_value
                nuc_poistion_values[nucIndex] = nuc_poistion_values[nucIndex] + adder
        
        #now record if the nuc position has a weghted bond
        for nucIndex in range(structure_list.nuc_count):
            is_weighted_bond=False
            if nuc_poistion_values[nucIndex] > struct_count:
                is_weighted_bond = True
            good_nucs_each_pos.append(is_weighted_bond)

        """
        for nucIndex in range(structure_list.nuc_count):
            nuc_value: float = float(nuc_poistion_values[nucIndex])
            #worked out this algotithm one night.. idk
            num_bonds_found:float = nuc_value / 2 + 1
            min_good_bonds: float = (num_bonds_found * 2) - ((nuc_value-num_bonds_found) * (-1))
            is_weighted_bond=False
            if num_bonds_found >= min_good_bonds and num_bonds_found > 0:
                is_weighted_bond = True
            good_nucs_each_pos.append(is_weighted_bond)
        """

        weighted_structure:str = ''
        for nucIndex in range(structure_list.nuc_count):
            is_bonded = good_nucs_each_pos[nucIndex]
            new_counter: collections.Counter = collections.Counter(nuc_pairs_comp_list[nucIndex])
            most_common_char: str= '.'
            if is_bonded is True:
                #most_common_char = '|'
                new_char:str = new_counter.most_common(2)[0][0]
                length = len(new_counter.most_common(2))
                if new_char == '.' and length > 1:
                    #then get second most common
                    new_char = new_counter.most_common(2)[1][0]
                most_common_char = new_char
            weighted_structure = weighted_structure + most_common_char

        weighted_structure: Sara2SecondaryStructure = Sara2SecondaryStructure(sequence=structure_list.sara_stuctures[0].sequence,
                                                                                structure=weighted_structure,)

        return weighted_structure