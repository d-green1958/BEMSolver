# Create a blade geometry class
from pandas import read_csv, errors
from numpy import pi as pi, unique, resize
from .aerofoil_lookup import AerofoilLookup


class BladeGeometry:
    def __init__(self):
        # known quantities 
        self.radial_distances = [] # {r}
        self.radial_differences = [] # {dr}
        self.chord_lengths = [] # c(r)
        self.twist_angles = [] # beta(r)
        self.drag_coeff = [] # C_l
        self.lift_coeff = [] # C_d
        self.aerofoil_type = [] # the aerofoil type at r
        self.B = 0 # number of blades
        
        # dictionary for different lookups for aerofoil data
        self.aerofoil_dict = {} # a dictionary: key - type name,
        # value aerofoil lookup
        
        # unknown quantities
        self.axial_inductance = [] # a
        self.tangential_inductance = [] # a'
        self.phi = [] # phi
        
        # calculated quantities
        self.angle_of_attack = [] # alpha
        self.solidity = []
        
    def set_number_of_blades(self, B: int):
        self.B = B
        print("Number of blades set.")        
        
        
    def update_solidity(self):
        if (self.B==0):
            raise Exception("B is 0. Please set the number of blades.")
        
        R = self.radial_distances[-1]
        temp = self.B/(2* pi * R)
        self.solidity = self.chord_lengths*temp/ self.radial_distances
        print("Solidity calculated.")
        
    
    def update_coeffs(self):
        for index, (r,c, alpha, type) in enumerate(zip(self.radial_distances, 
                                        self.chord_lengths,
                                        self.angle_of_attack,
                                        self.aerofoil_type)):
           
            self.lift_coeff[index], self.drag_coeff[index] = (self.aerofoil_dict.get(type)).get_coeff(alpha)

            
            
    
    def read_configuration(self, file_path: str):
        config_name = (file_path.split('/')[-1]).split(".case")[0]
        try:
            config = read_csv(file_path, delim_whitespace=True,
                        comment="!", names=['r', 'beta', 'dr', 'chord', 'type'])
            print("Configuration read succesfully.")
        except errors.EmptyDataError:
            print("Empty file:", self.aerofoil_name)
        except FileNotFoundError:
            print(f"File not found:", self.aerofoil_name)
            
        # add the data from the configuration files
        self.radial_distances = config.r
        self.radial_differences = config.dr
        self.twist_angles = config.beta
        self.chord_lengths = config.chord
        self.aerofoil_type = config.type
        
        # resize some stuff
        length = len(self.radial_distances)
        self.lift_coeff = [0] * length
        self.drag_coeff = [0] * length
        self.angle_of_attack = [0] * length
                
        # now add the lookups for the aerofoil data
        for aerofoil_type in unique(config.type):
            print("Adding: " + aerofoil_type)
            self.aerofoil_dict[aerofoil_type] = AerofoilLookup(aerofoil_type)
            
        print("Configuration added: ", config_name)
        print("Number of sections:  ", length)
        print("Blade Length:  ", max(self.radial_distances))
    
    def calculate_angle_of_attack(self):
        self.angle_of_attack = self.phi - self.twist_angles
        
    def update_factors(self, axial_inductance, tangential_inductance):
        self.axial_inductance = axial_inductance
        self.tangential_inductance = tangential_inductance
        print("Inductance factors updated.")
        
    def update_phi(self, phi):
        self.phi = phi
        print("Phi updated.")
    
    
