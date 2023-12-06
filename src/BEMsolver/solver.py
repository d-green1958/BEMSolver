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
    
    
    
    
    
if __name__ == "__main__":
    print("DONT PRINT")
    # blade = BladeGeometry()
