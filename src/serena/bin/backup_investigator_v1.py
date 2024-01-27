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


from serena.utilities.ensemble_structures import Sara2SecondaryStructure
from serena.utilities.ensemble_structures import Sara2StructureList
from serena.analysis.investigator import ComparisonEvalResults
from serena.analysis.investigator import RatioResults
from serena.utilities.comparison_structures import ComparisonNucResults
from serena.utilities.comparison_structures import ComparisonNucCounts
from serena.utilities.local_minima_variation import ComparisonLMVResponse
from serena.utilities.local_minima_variation import ComparisonLMV
from serena.analysis.investigator import LMVAssertionResult
from serena.utilities.ensemble_variation import EV
from serena.analysis.scoring import BasicScoreResults
from serena.analysis.scoring import AdvancedScoreResults
from serena.interfaces.Sara2_API_Python3 import DesignInformation
from serena.interfaces.Sara2_API_Python3 import WetlabData

class Nut_Attributes(Enum):
	Investigator = "investigator_db"
	Scores = "scores_db"
	DesignParameters = "design_info_db"


class PNASData(Nut):

	def __init__(self, working_folder:Path, var_name:str, use_db:bool = False) -> None:
		super().__init__(enum_list=Nut_Attributes,
			use_db=True,
			db=None,
			var_name=var_name,
			working_folder=working_folder)


		self.investigator_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="comparison_eval_result_db",
			atr_type=['ComparisonEvalResults', 'RatioResults']))

		self.investigator_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="comp_nuc_counts_db",
			atr_type=['ComparisonNucResults', 'ComparisonNucCounts']))

		self.investigator_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="lmv_values_db",
			atr_type=['ComparisonLMVResponse', 'ComparisonLMV', 'EV']))

		self.investigator_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="lmv_assertions_db",
			atr_type=['LMVAssertionResult']))

		self.investigator_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="num_groups_db",
			atr_type=int))

		self.investigator_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="total_structures_ensemble_db",
			atr_type=int))

		self.scores_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="basic_scores_db",
			atr_type=['BasicScoreResults']))

		self.scores_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="advanced_scores_db",
			atr_type=['AdvancedScoreResults']))

		self.scores_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="number_structures_db",
			atr_type=int))

		self.design_info_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="design_info_db",
			atr_type=['DesignInformation']))

		self.design_info_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="wetlab_data_db",
			atr_type=['WetlabData']))

class DesignParameters(CustomAttribute):
	def __init__(self, parent: Any, current:Any, save_value:bool) -> None:
		self.parent = parent
		self.current = current
		self.do_save = save_value

	@property
	def design_info(self)->DesignInformation:
		self.parent.nut_filter.yaml_operations.yaml.register_class(DesignInformation)
		return self.parent.design_info_db

	@design_info.setter
	def design_info(self, value:DesignInformation):
		if isinstance(value, DesignInformation) == False:
			raise ValueError("Invalid value assignment")
		self.parent.nut_filter.yaml_operations.yaml.register_class(DesignInformation)
		self.parent.design_info_db = value


	@property
	def wetlab_data(self)->WetlabData:
		self.parent.nut_filter.yaml_operations.yaml.register_class(WetlabData)
		return self.parent.wetlab_data_db

	@wetlab_data.setter
	def wetlab_data(self, value:WetlabData):
		if isinstance(value, WetlabData) == False:
			raise ValueError("Invalid value assignment")
		self.parent.nut_filter.yaml_operations.yaml.register_class(WetlabData)
		self.parent.wetlab_data_db = value


class Investigator(CustomAttribute):
	def __init__(self, parent: Any, current:Any, save_value:bool) -> None:
		self.parent = parent
		self.current = current
		self.do_save = save_value

	@property
	def comparison_eval_result(self)->ComparisonEvalResults:
		self.parent.nut_filter.yaml_operations.yaml.register_class(ComparisonEvalResults)
		self.parent.nut_filter.yaml_operations.yaml.register_class(RatioResults)
		return self.parent.comparison_eval_result_db

	@comparison_eval_result.setter
	def comparison_eval_result(self, value:ComparisonEvalResults):
		if isinstance(value, ComparisonEvalResults) == False:
			raise ValueError("Invalid value assignment")
		self.parent.nut_filter.yaml_operations.yaml.register_class(ComparisonEvalResults)
		self.parent.nut_filter.yaml_operations.yaml.register_class(RatioResults)
		self.parent.comparison_eval_result_db = value


	@property
	def comp_nuc_counts(self)->ComparisonNucResults:
		self.parent.nut_filter.yaml_operations.yaml.register_class(ComparisonNucResults)
		self.parent.nut_filter.yaml_operations.yaml.register_class(ComparisonNucCounts)
		return self.parent.comp_nuc_counts_db

	@comp_nuc_counts.setter
	def comp_nuc_counts(self, value:ComparisonNucResults):
		if isinstance(value, ComparisonNucResults) == False:
			raise ValueError("Invalid value assignment")
		self.parent.nut_filter.yaml_operations.yaml.register_class(ComparisonNucResults)
		self.parent.nut_filter.yaml_operations.yaml.register_class(ComparisonNucCounts)
		self.parent.comp_nuc_counts_db = value


	@property
	def lmv_values(self)->ComparisonLMVResponse:
		self.parent.nut_filter.yaml_operations.yaml.register_class(ComparisonLMVResponse)
		self.parent.nut_filter.yaml_operations.yaml.register_class(ComparisonLMV)
		self.parent.nut_filter.yaml_operations.yaml.register_class(EV)
		return self.parent.lmv_values_db

	@lmv_values.setter
	def lmv_values(self, value:ComparisonLMVResponse):
		if isinstance(value, ComparisonLMVResponse) == False:
			raise ValueError("Invalid value assignment")
		self.parent.nut_filter.yaml_operations.yaml.register_class(ComparisonLMVResponse)
		self.parent.nut_filter.yaml_operations.yaml.register_class(ComparisonLMV)
		self.parent.nut_filter.yaml_operations.yaml.register_class(EV)
		self.parent.lmv_values_db = value


	@property
	def lmv_assertions(self)->LMVAssertionResult:
		self.parent.nut_filter.yaml_operations.yaml.register_class(LMVAssertionResult)
		return self.parent.lmv_assertions_db

	@lmv_assertions.setter
	def lmv_assertions(self, value:LMVAssertionResult):
		if isinstance(value, LMVAssertionResult) == False:
			raise ValueError("Invalid value assignment")
		self.parent.nut_filter.yaml_operations.yaml.register_class(LMVAssertionResult)
		self.parent.lmv_assertions_db = value


	@property
	def num_groups(self)->int:
		return self.parent.num_groups_db

	@num_groups.setter
	def num_groups(self, value:int):
		if isinstance(value, int) == False:
			raise ValueError("Invalid value assignment")
		self.parent.num_groups_db = value


	@property
	def total_structures_ensemble(self)->int:
		return self.parent.total_structures_ensemble_db

	@total_structures_ensemble.setter
	def total_structures_ensemble(self, value:int):
		if isinstance(value, int) == False:
			raise ValueError("Invalid value assignment")
		self.parent.total_structures_ensemble_db = value


class Scores(CustomAttribute):
	def __init__(self, parent: Any, current:Any, save_value:bool) -> None:
		self.parent = parent
		self.current = current
		self.do_save = save_value

	@property
	def basic_scores(self)->BasicScoreResults:
		self.parent.nut_filter.yaml_operations.yaml.register_class(BasicScoreResults)
		return self.parent.basic_scores_db

	@basic_scores.setter
	def basic_scores(self, value:BasicScoreResults):
		if isinstance(value, BasicScoreResults) == False:
			raise ValueError("Invalid value assignment")
		self.parent.nut_filter.yaml_operations.yaml.register_class(BasicScoreResults)
		self.parent.basic_scores_db = value


	@property
	def advanced_scores(self)->AdvancedScoreResults:
		self.parent.nut_filter.yaml_operations.yaml.register_class(AdvancedScoreResults)
		return self.parent.advanced_scores_db

	@advanced_scores.setter
	def advanced_scores(self, value:AdvancedScoreResults):
		if isinstance(value, AdvancedScoreResults) == False:
			raise ValueError("Invalid value assignment")
		self.parent.nut_filter.yaml_operations.yaml.register_class(AdvancedScoreResults)
		self.parent.advanced_scores_db = value


	@property
	def number_structures(self)->int:
		return self.parent.number_structures_db

	@number_structures.setter
	def number_structures(self, value:int):
		if isinstance(value, int) == False:
			raise ValueError("Invalid value assignment")
		self.parent.number_structures_db = value


class ArchiveInvestigator(PNASData):

	def __init__(self, working_folder:str, var_name:str, use_db:bool = False) -> None:
		super().__init__(use_db=use_db,
			var_name=var_name,
			working_folder=Path(working_folder))


		self._investigator: Investigator = Investigator(save_value=True,
			current=None,
			parent=self.investigator_db)

		self._scores: Scores = Scores(save_value=True,
			current=None,
			parent=self.scores_db)

		self._design_info: DesignParameters = DesignParameters(save_value=True,
			current=None,
			parent=self.design_info_db)

	@property
	def investigator(self)->Investigator:
		return self._investigator

	@investigator.setter
	def investigator(self, struct:Investigator):
		if isinstance(struct, Investigator) == False:
			raise ValueError("Invalid value assignment")
		self._investigator = struct


	@property
	def scores(self)->Scores:
		return self._scores

	@scores.setter
	def scores(self, struct:Scores):
		if isinstance(struct, Scores) == False:
			raise ValueError("Invalid value assignment")
		self._scores = struct


	@property
	def design_info(self)->DesignParameters:
		return self._design_info

	@design_info.setter
	def design_info(self, struct:DesignParameters):
		if isinstance(struct, DesignParameters) == False:
			raise ValueError("Invalid value assignment")
		self._design_info = struct


