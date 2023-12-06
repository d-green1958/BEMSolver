# Create a blade geometry class

from numpy import pi as pi


class BladeGeometry:
    def __init__(self):
        # known quantities 
        self.radial_distances = [] # {r}
        self.chord_lengths = [] # c(r)
        self.drag_values = [] # C_d(r)
        self.lift_values = [] # C_l(r)
        self.twist_angles = [] # beta(r)
        self.B = 0 # number of blades
        
        # unknown quantities
        self.axial_inductance = [] # a
        self.tangential_inductance = [] # a'
        self.phi = [] # phi
        
        # calculated quantities
        self.angle_of_attack = [] # alpha
        self.solidity = []
        
        # look up table for the drag and chord lengths
        self.aerofoil_lookup_table = []
        
    def update_factors(self, axial_inductance, tangential_inductance):
        self.axial_inductance = axial_inductance
        self.tangential_inductance = tangential_inductance
        print("Inductance factors updated.")
        
    def update_phi(self, phi):
        self.phi = phi
        print("Phi updated.")
    
    def append_data_point(self, radius, chord_length, drag, lift, twist_angle):
        self.radial_distances.append(radius)
        self.chord_lengths.append(chord_length)
        self.drag_values.append(drag)
        self.lift_values.append(lift)
        self.twist_angles.append(twist_angle)
        print("Data point added.")
        
    def add_dataset(self, chord_length, drag, lift, twist_angle):
        self.chord_lengths = chord_length
        self.drag_values = drag
        self.lift_values = lift
        self.twist_angles = twist_angle
        print("Data set added.")
    
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
    
    def display_data(self):
        print("Radial Distance\tChord Length\tDrag\tLift\tTwist Angle")
        for i in range(len(self.radial_distances)):
            print(f"{self.radial_distances[i]}\t\t{self.chord_lengths[i]}\t\t{self.drag_values[i]}\t{self.lift_values[i]}\t{self.twist_angles[i]}")
    
