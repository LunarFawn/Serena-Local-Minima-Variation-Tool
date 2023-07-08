"""
Stuff I spun up just for eternacon presentation
"""


from dataclasses import dataclass
from serena.utilities.Sara2_API_Python3 import Sara2API

@dataclass
class PNASWetLab():
    pass

class run_eternacon():

    def __init__(self) -> None:
        pass

    def run(self):
        """
        run this
        """
        
        pnas_path:str = ''
        pnas_round101_sheet:str = ''
        
        new_sara:Sara2API = Sara2API()

        new_sara.ProcessLab()

