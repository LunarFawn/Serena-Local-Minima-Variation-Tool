"""
File that defines the main RNA sequence data
"""


from enum import Enum
from attrs import define, field
from collections import namedtuple
from typing import List, Dict, Any,TypeVar, Type
from pathlib import Path

from data_squirrel.config.dynamic_data_nut import (
	Nut,
	Value,
	GenericAttribute,
	AtrClass,
	CustomAttribute
)


class Nut_Attributes(Enum):
	Sara2StructureList = "structs_db"
	Sara2SecondaryStructure = "sara2_db"


class RNAStrand(Nut):

	def __init__(self, working_folder:Path, var_name:str, use_db:bool = False) -> None:
		super().__init__(enum_list=Nut_Attributes,
			use_db=True,
			db=None,
			var_name=var_name,
			working_folder=working_folder)


		self.structs_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="mfe_structure_db",
			atr_type=str))

		self.structs_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="mfe_free_energy_db",
			atr_type=float))

		self.structs_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="mfe_stack_energy_db",
			atr_type=float))

		self.structs_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="nuc_count_db",
			atr_type=int))

		self.structs_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="sara_stuctures_db",
			atr_type=['Sara2SecondaryStructure', 'CLASS']))

		self.structs_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="max_free_energy_db",
			atr_type=float))

		self.structs_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="min_free_energy_db",
			atr_type=float))

		self.structs_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="max_stack_energy_db",
			atr_type=float))

		self.structs_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="min_stack_energy_db",
			atr_type=float))

		self.structs_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="num_structures_db",
			atr_type=int))

		self.structs_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="free_energy_span_db",
			atr_type=float))

		self.structs_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="stack_energy_span_db",
			atr_type=float))

		self.structs_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="weighted_structure_db",
			atr_type=str))

		self.sara2_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="sequence_db",
			atr_type=str))

		self.sara2_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="structure_db",
			atr_type=str))

		self.sara2_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="free_energy_db",
			atr_type=float))

		self.sara2_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="stack_energy_db",
			atr_type=float))

class Sara2SecondaryStructure(CustomAttribute):
	def __init__(self, parent: Any, current:Any, save_value:bool) -> None:
		self.parent = parent
		self.current = current
		self.do_save = save_value

	@property
	def sequence(self)->str:
		return self.parent.sequence_db

	@sequence.setter
	def sequence(self, value:str):
		if isinstance(value, str) == False:
			raise ValueError("Invalid value assignment")
		self.parent.sequence_db = value


	@property
	def structure(self)->str:
		return self.parent.structure_db

	@structure.setter
	def structure(self, value:str):
		if isinstance(value, str) == False:
			raise ValueError("Invalid value assignment")
		self.parent.structure_db = value


	@property
	def free_energy(self)->float:
		return self.parent.free_energy_db

	@free_energy.setter
	def free_energy(self, value:float):
		if isinstance(value, float) == False:
			raise ValueError("Invalid value assignment")
		self.parent.free_energy_db = value


	@property
	def stack_energy(self)->float:
		return self.parent.stack_energy_db

	@stack_energy.setter
	def stack_energy(self, value:float):
		if isinstance(value, float) == False:
			raise ValueError("Invalid value assignment")
		self.parent.stack_energy_db = value


class Sara2StructureList(CustomAttribute):
	def __init__(self, parent: Any, current:Any, save_value:bool) -> None:
		self.parent = parent
		self.current = current
		self.do_save = save_value

	@property
	def mfe_structure(self)->str:
		return self.parent.mfe_structure_db

	@mfe_structure.setter
	def mfe_structure(self, value:str):
		if isinstance(value, str) == False:
			raise ValueError("Invalid value assignment")
		self.parent.mfe_structure_db = value


	@property
	def mfe_free_energy(self)->float:
		return self.parent.mfe_free_energy_db

	@mfe_free_energy.setter
	def mfe_free_energy(self, value:float):
		if isinstance(value, float) == False:
			raise ValueError("Invalid value assignment")
		self.parent.mfe_free_energy_db = value


	@property
	def mfe_stack_energy(self)->float:
		return self.parent.mfe_stack_energy_db

	@mfe_stack_energy.setter
	def mfe_stack_energy(self, value:float):
		if isinstance(value, float) == False:
			raise ValueError("Invalid value assignment")
		self.parent.mfe_stack_energy_db = value


	@property
	def nuc_count(self)->int:
		return self.parent.nuc_count_db

	@nuc_count.setter
	def nuc_count(self, value:int):
		if isinstance(value, int) == False:
			raise ValueError("Invalid value assignment")
		self.parent.nuc_count_db = value


	@property
	def sara_stuctures(self)->List[Sara2SecondaryStructure]:
		new_list:List[Sara2SecondaryStructure] = []

		temp_list = self.parent.sara_stuctures_db

		for index, item in enumerate(temp_list):
			temp_parent:str = f'entry{index}'
			new_parent:CustomAttribute = self.parent.new_attr(GenericAttribute(atr_class=AtrClass.PARENT,
				attribute=temp_parent,
				atr_type=None))

			temp_attr2 = getattr(new_parent, temp_parent)
			temp_name:str = 'temp'
			setattr(self, temp_name, Sara2SecondaryStructure(parent=temp_attr2,
				current=None,
				save_value=True))

			for key, value in item.items():
				setattr(temp_attr2, key, value)
			new_primary_struct:Sara2SecondaryStructure = getattr(self, temp_name)
			new_list.append(new_primary_struct)

		return new_list

	@sara_stuctures.setter
	def sara_stuctures(self, value:List[Sara2SecondaryStructure]):
		if isinstance(value, list) == False:
			raise ValueError("Invalid value assignment")
		if len(value) < 1:
			raise Exception("Empty lists not allowed")

		for item in value:
			if isinstance(item, Sara2SecondaryStructure) == False:
				raise ValueError("Invalid value assignment")
		self.parent.sara_stuctures_db = value


	@property
	def max_free_energy(self)->float:
		return self.parent.max_free_energy_db

	@max_free_energy.setter
	def max_free_energy(self, value:float):
		if isinstance(value, float) == False:
			raise ValueError("Invalid value assignment")
		self.parent.max_free_energy_db = value


	@property
	def min_free_energy(self)->float:
		return self.parent.min_free_energy_db

	@min_free_energy.setter
	def min_free_energy(self, value:float):
		if isinstance(value, float) == False:
			raise ValueError("Invalid value assignment")
		self.parent.min_free_energy_db = value


	@property
	def max_stack_energy(self)->float:
		return self.parent.max_stack_energy_db

	@max_stack_energy.setter
	def max_stack_energy(self, value:float):
		if isinstance(value, float) == False:
			raise ValueError("Invalid value assignment")
		self.parent.max_stack_energy_db = value


	@property
	def min_stack_energy(self)->float:
		return self.parent.min_stack_energy_db

	@min_stack_energy.setter
	def min_stack_energy(self, value:float):
		if isinstance(value, float) == False:
			raise ValueError("Invalid value assignment")
		self.parent.min_stack_energy_db = value


	@property
	def num_structures(self)->int:
		return self.parent.num_structures_db

	@num_structures.setter
	def num_structures(self, value:int):
		if isinstance(value, int) == False:
			raise ValueError("Invalid value assignment")
		self.parent.num_structures_db = value


	@property
	def free_energy_span(self)->float:
		return self.parent.free_energy_span_db

	@free_energy_span.setter
	def free_energy_span(self, value:float):
		if isinstance(value, float) == False:
			raise ValueError("Invalid value assignment")
		self.parent.free_energy_span_db = value


	@property
	def stack_energy_span(self)->float:
		return self.parent.stack_energy_span_db

	@stack_energy_span.setter
	def stack_energy_span(self, value:float):
		if isinstance(value, float) == False:
			raise ValueError("Invalid value assignment")
		self.parent.stack_energy_span_db = value


	@property
	def weighted_structure(self)->str:
		return self.parent.weighted_structure_db

	@weighted_structure.setter
	def weighted_structure(self, value:str):
		if isinstance(value, str) == False:
			raise ValueError("Invalid value assignment")
		self.parent.weighted_structure_db = value


class BackupSerena(RNAStrand):

	def __init__(self, working_folder:str, var_name:str, use_db:bool = False) -> None:
		super().__init__(use_db=use_db,
			var_name=var_name,
			working_folder=Path(working_folder))


		self._structs: Sara2StructureList = Sara2StructureList(save_value=True,
			current=None,
			parent=self.structs_db)

		self._sara2: Sara2SecondaryStructure = Sara2SecondaryStructure(save_value=True,
			current=None,
			parent=self.sara2_db)

	@property
	def structs(self)->Sara2StructureList:
		return self._structs

	@structs.setter
	def structs(self, struct:Sara2StructureList):
		if isinstance(struct, Sara2StructureList) == False:
			raise ValueError("Invalid value assignment")
		self._structs = struct


	@property
	def sara2(self)->Sara2SecondaryStructure:
		return self._sara2

	@sara2.setter
	def sara2(self, struct:Sara2SecondaryStructure):
		if isinstance(struct, Sara2SecondaryStructure) == False:
			raise ValueError("Invalid value assignment")
		self._sara2 = struct


