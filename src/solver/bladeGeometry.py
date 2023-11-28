# Create a blade geometry class

class BladeGeometry:
    def __init__(self, radial_distances):
        # known quantities 
        self.radial_distances = radial_distances
        self.chord_lengths = []
        self.drag_values = []
        self.lift_values = []
        self.twist_angles = []
        
        # unknown quantities
        self.axial_inductance = []
        self.tangential_inductance = []
        self.phi = []
        
        # calculated quantities
        self.angle_of_attack
        
    def update_factors(self, axial_inductance, tangential_inductance):
        self.axial_inductance = axial_inductance
        self.tangential_inductance = tangential_inductance
        
    def update_phi(self, phi):
        self.phi = phi
    
    def append_data_point(self, chord_length, drag, lift, twist_angle):
        self.chord_lengths.append(chord_length)
        self.drag_values.append(drag)
        self.lift_values.append(lift)
        self.twist_angles.append(twist_angle)
        
    def add_data(self, chord_length, drag, lift, twist_angle):
        self.chord_lengths = chord_length
        self.drag_values = drag
        self.lift_values = lift
        self.twist_angles = twist_angle
        
    def display_data(self):
        print("Radial Distance\tChord Length\tDrag\tLift\tTwist Angle")
        for i in range(len(self.radial_distances)):
            print(f"{self.radial_distances[i]}\t\t{self.chord_lengths[i]}\t\t{self.drag_values[i]}\t{self.lift_values[i]}\t{self.twist_angles[i]}")
    
