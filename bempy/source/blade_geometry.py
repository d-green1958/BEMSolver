# Create a blade geometry class
from pandas import read_csv, errors
from numpy import pi as pi, unique, deg2rad, rad2deg
from .aerofoil_lookup import AerofoilLookup



class BladeGeometry:
    def __init__(self, silent_mode = False):
        # known quantities 
        self.number_of_nodes = 0 # number of blade sections
        self.radial_distances = [] # {r}
        self.radial_differences = [] # {dr}
        self.chord_lengths = [] # c(r)
        self.twist_angles = [] # beta(r)
        self.drag_coeff = [] # C_l
        self.lift_coeff = [] # C_d
        self.alpha_0 = [] # AOA for no lift
        self.aerofoil_type = [] # the aerofoil type at r
        self.B = 3 # number of blades
        
        # dictionary for different lookups for aerofoil data
        self.aerofoil_dict = {} # a dictionary: key - type name,
        # value aerofoil lookup
        
        # unknown quantities
        self.axial_inductance = [] # a
        self.tangential_inductance = [] # a'
        
        self.axial_inductance_prev = [] # a
        self.tangential_inductance_prev = [] # a'
        
        self.phi = [] # phi
        self.R = 0 # blade tip speed radius
        self.R_hub = 0 # the blade hub radius
        self.blade_length = 0 # blade length
        
        # calculated quantities
        self.angle_of_attack = [] # alpha
        self.solidity = []    
        self.chord_solidity = []
        
        self.silent_mode = silent_mode
    
    def print_details(self):
        print(f"{'Node':<10} {'r [m]':<10} {'dr [m]':<10} {'c(r) [m]':<10} {'beta(r) [rad]':<15} {'type name':<15}")
        for node in range(self.number_of_nodes):
            print(f"{node:<10.0f} {self.radial_distances[node]:<10.4f} {self.radial_differences[node]:<10.4f} {self.chord_lengths[node]:<10.4f} {self.twist_angles[node]:<15.4f} {self.aerofoil_type[node]:<15}")
        
    def calculate_solidity(self):
        if (self.B==0):
            raise Exception("B is 0. Please set the number of blades.")
        
        temp = self.B/(2* pi * self.R)
        self.solidity = self.chord_lengths*temp/ self.radial_distances
        if not self.silent_mode:
            print("Solidity calculated (B = " + str(self.B) + ").")
        
    def calculate_chord_solidity(self):
        temp = self.B/(2* pi)
        self.chord_solidity = temp*self.chord_lengths/self.radial_distances
        
    def read_configuration(self, file_path: str):
        config_name = (file_path.split('/')[-1]).split(".case")[0]
        try:
            config = read_csv(file_path, delim_whitespace=True,
                        comment="!", names=['r', 'beta', 'dr', 'chord', 'type'])
            if not self.silent_mode:
                print("Configuration read succesfully.")
        except errors.EmptyDataError:
            print("Empty file:", self.aerofoil_name)
        except FileNotFoundError():
            print(f"File not found:", self.aerofoil_name)
            
        # add the data from the configuration files
        self.radial_distances = config.r
        self.radial_differences = config.dr
        self.twist_angles = deg2rad(config.beta)
        self.chord_lengths = config.chord
        self.aerofoil_type = config.type
        
        # resize some stuff
        length = len(self.radial_distances)
        self.number_of_nodes = length
        self.lift_coeff = [0] * length
        self.drag_coeff = [0] * length
        self.angle_of_attack = [0] * length
        
        self.axial_inductance = [0] * length
        self.tangential_inductance = [0] * length
        self.axial_inductance_prev = [0] * length
        self.tangential_inductance_prev = [0] * length
        self.phi = [0] * length
        self.alpha_0 = [0] * length
                
        # now add the lookups for the aerofoil data
        for aerofoil_type in unique(config.type):
            if not self.silent_mode:
                print("Adding: " + aerofoil_type)
            self.aerofoil_dict[aerofoil_type] = AerofoilLookup(aerofoil_type)
    
        self.R = list(self.radial_distances)[-1] + 0.5*list(self.radial_differences)[-1]
        self.R_hub = list(self.radial_distances)[0] - 0.5*list(self.radial_differences)[0]
        self.blade_length = round(sum(self.radial_differences),3)
        
        for node in range(self.number_of_nodes):
            self.alpha_0[node] = (self.aerofoil_dict[self.aerofoil_type[node]]).alpha_0
            
        if not self.silent_mode:
            print("Configuration added: ", config_name)
            print("Number of sections:  ", length)
            print("Blade Length:  ", self.blade_length)
            print("Blade tip radius:  ", self.R)
            print()
            self.print_details()
    
    def update_angle_of_attack(self, node = -1):
        if node == -1:
            self.angle_of_attack = rad2deg(self.phi - self.twist_angles)
        else:
            self.angle_of_attack[node] = rad2deg(self.phi[node] - self.twist_angles[node])
        
    def update_drag_and_lift(self, node = -1):
        if node == -1:
            for ind in range(self.number_of_nodes):
                self.lift_coeff[ind], self.drag_coeff[ind] = (self.aerofoil_dict[self.aerofoil_type[ind]]).get_coeff(self.angle_of_attack[ind])
        else:
            self.lift_coeff[node], self.drag_coeff[node] = (self.aerofoil_dict[self.aerofoil_type[node]]).get_coeff(self.angle_of_attack[node])

            
            
    
    
