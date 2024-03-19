from pandas import read_csv, errors
from os import path
from numpy import interp



# class used to lookup and import aerofoil data
class AerofoilLookup:
    def __init__(self, aerofoil_name: str, dir: str = "/../aerofoil_data/"):
        # create empty arrays
        self.angle_of_attack = []
        self.C_drag = []
        self.C_lift = []
        
        
        self.drag_dict = {}
        self.lift_dict = {}
        
        # path to the aerofoil data directory
        dir_path = path.dirname(__file__)
        self.path = dir_path + dir
        
        # path to specific file
        self.aerofoil_name = aerofoil_name
        
        # load the data
        self.load_data()
    
    def __str__(self):
        return "AerofoilLookup: data - " + self.aerofoil_name
        
        
    def load_data(self):
        try:
            data = read_csv(self.path + self.aerofoil_name +".dat", delim_whitespace=True,
                        comment="!", names=['Alpha', 'Cl', 'Cd', 'Cm'])
        except errors.EmptyDataError:
            print("Empty file:", self.aerofoil_name)
        except FileNotFoundError:
            print(f"File not found:", self.aerofoil_name)
            
        self.angle_of_attack = data.Alpha
        
        # Again legacy
        # self.lift_dict = dict(zip(data.Alpha, data.Cl))
        # self.drag_dict = dict(zip(data.Alpha, data.Cd))
        
        self.C_lift = data.Cl
        self.C_drag = data.Cd
    
    def round_alpha(self, test_value):
        return min(self.angle_of_attack, key = lambda x: abs(x - test_value))
    
    def get_coeff(self, alpha: float) -> list[float]:
        
        # legacy code when alpha was rounded to nearest value
        # used to find the nearest value of alpha in the list
        # alpha = self.round_alpha(alpha)
        # return [self.lift_dict[alpha], self.drag_dict[alpha]]
        
        # instead use linear interpolation between alpha values
        
        Cl = interp(alpha, self.angle_of_attack, self.C_lift)
        Cd = interp(alpha, self.angle_of_attack, self.C_drag)
        
        return [Cl, Cd]
        
        
        
        
        

    
    