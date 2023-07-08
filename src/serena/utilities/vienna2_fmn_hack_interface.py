
import subprocess
from dataclasses import dataclass

@dataclass
class RNAFoldResponse():
    sequence:str
    structure:str
    energy:float

class Vienna2FMNInterface():
    """
    class to interface with vienna2_fmn hack by Elnando888
    """
    
    def __init__(self) -> None:
        pass

    def rnafold_fmn(self,input_sequence:str, fmn_amount:int = 200)->RNAFoldResponse:
        command= ["RNAfold", "--ligand", f'FMN:{fmn_amount}']
        sequence: str = ''
        structure:str = ''
        energy:float = 0
        process = subprocess.run(command, input=input_sequence, encoding="utf-8", capture_output=True)
        raw_result = process.stdout.split('\n')
        if len(raw_result) > 0:
            sequence = raw_result[0]
            raw_struct = raw_result[1].split(' ')
            if len(raw_struct) > 0:
                structure = raw_struct[0]
                energy_str:str = raw_struct[1]
                energy_str = energy_str.strip('(')
                energy_str = energy_str.strip(')')
                try:
                    energy = float(energy_str)
                except ValueError:
                    raise Exception("not a number")
        
        response: RNAFoldResponse = RNAFoldResponse(sequence=sequence,
                                                    structure=structure,
                                                    energy=energy)
        return response
                        
    def run_via_ssh():
        pass



new_viena: Vienna2FMNInterface = Vienna2FMNInterface()
seq = 'GCCAUCGCAUGAGGAUAUGCUCGGGUUUCCCGAGCAGAAGGCAUGUCACAAGACAUGAGGAUCACCCAUGUAGUUAAGAUGGCA'
result = new_viena.rnafold_fmn(seq)
print(result)