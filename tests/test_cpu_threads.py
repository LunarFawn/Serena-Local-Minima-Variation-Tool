import pytest


from ..src.serena.utilities.cpu_thread_manager import CPU_Data
from ..src.serena.utilities.data import SwitchAnalysisDataToProcess

@pytest.fixture
def empty_cpu_data():
    return CPU_Data()

@pytest.fixture
def filled_switch_data():    
    return SwitchAnalysisDataToProcess(temperature=10,
                                        fmn_struct = '...',
                                        fmn_struct_free_energy=20,
                                        sequence = '(((',
                                        span = 30,
                                        units = 40)

def test_cpu_data_default(empty_cpu_data:CPU_Data):
    assert(empty_cpu_data.temperature_data.sequence=='')
    assert(empty_cpu_data.temperature_data.span==-1)
    assert(empty_cpu_data.temperature_data.units==-1)
    assert(empty_cpu_data.temperature_data.fmn_struct=='')
    assert(empty_cpu_data.temperature_data.fmn_struct_free_energy==-1)
    assert(empty_cpu_data.temperature_data.temperature==-1)


def test_fill_cpu_data(filled_switch_data:SwitchAnalysisDataToProcess):
    assert(filled_switch_data.temperature==10)
    assert(filled_switch_data.fmn_struct == '...')
    assert(filled_switch_data.fmn_struct_free_energy==20)
    assert(filled_switch_data.sequence == '(((')
    assert(filled_switch_data.span == 30)
    assert(filled_switch_data.units == 40)

def test_replace_cpu_data(empty_cpu_data:CPU_Data, filled_switch_data:SwitchAnalysisDataToProcess):
    empty_cpu_data.temperature_data = filled_switch_data
    assert(empty_cpu_data.temperature_data.temperature==10)
    assert(empty_cpu_data.temperature_data.fmn_struct == '...')
    assert(empty_cpu_data.temperature_data.fmn_struct_free_energy==20)
    assert(empty_cpu_data.temperature_data.sequence == '(((')
    assert(empty_cpu_data.temperature_data.span == 30)
    assert(empty_cpu_data.temperature_data.units == 40)



