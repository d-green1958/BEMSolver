from .blade_geometry import BladeGeometry
from math import atan 

def calculate_phi_array(instance: BladeGeometry, inflow_velocity: float,
                        rotational_velocity: float) -> list[float]:
    return atan(instance.radial_distances * rotational_velocity *
                    (1 + instance.tangential_inductance) / 
                    (inflow_velocity * (1 - instance.axial_inductance)))
    
def calculate_axial_inductance(instance: BladeGeometry) -> list[float]:
    # finish this
    print("finish this")
    
    
class Problem:
    def __init__(self):
        self.blade = BladeGeometry()
    
    def add_configuration(self, config_path: str):
        boarder = "#"*10
        print(boarder + "  CONFIG  " + boarder)
        self.blade.read_configuration(config_path)
        print()
        print()
        
    def apply_ICs(self, phi_IC: list[float]):
        self.blade.phi = phi_IC
        self.blade.calculate_angle_of_attack()
        self.blade.update_coeffs()
    
    def run(self, tol: float, iter_max: int):
        boarder = "#"*10
        print(boarder + "  SOLVE  " + boarder)

        iter = 0
        print()
        while iter<iter_max:
            # finish iteration m ethod
            
            # need some
            # if err < tol:
            #     break
            iter += 1
        
        print()
        print()
    