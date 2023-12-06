from pandas import read_csv
from os import path


# class used to lookup and import aerofoil data
class AerofoilLookup:
    def __init__(self, aerofoil_name: str, dir: str = "/aerofoil_data/"):
        self.angle_of_attack = []
        self.C_drag = []
        self.C_lift = []
        
        dir_path = path.dirname(__file__)
        
        self.path = dir_path + dir
        self.aerofoil_name = aerofoil_name
        
        
    
    def change_dir(self, dir: str):
        self.path = dir
    
    def change_name(self, file_name: str):
        self.file_name = file_name 
        
    def get_data(self) -> list[float]:
        data = read_csv(self.path + self.aerofoil_name +".dat", delim_whitespace=True,
                        comment="!", names=['Alpha', 'Cl', 'Cd', 'Cm'])
        print("opened:", self.path + self.aerofoil_name +".dat")
        return [data.Alpha, data.Cl, data.Cd]

    
    