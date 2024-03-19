

# class used to present solution results in an organised manner
class Result:
    def __init__(self) -> None:
        
        # simulation parameters
        self.wind_speed = None
        self.rot_speed = None
        self.tip_speed_ratio = None
        self.rotor_radius = None
        
        # blade data
        self.sectional_positions = None
        self.sectional_chord_lengths = None
        self.sectional_lengths = None
        
        # sectional results
        self.sectional_torque = None
        self.sectional_thrust = None
        
        # integrated results
        self.power = None
        self.thrust = None
        self.torque = None
        
        # inductance factors and inflow angle
        self.sectional_axial_inductance_factor = None
        self.sectional_tangential_inductance_factor = None
        self.sectional_phi = None
        
        # simulation performance
        self.sectional_convergence = None
        self.sectional_iterations = None
        self.sectional_final_error = None
        
        self.power_coefficient = None
        self.torque_coefficient = None
        self.thrust_coefficient = None
        
        # possibility to add time varying components i.e. how axial inductance factors change with iterations and errors and other things like this
        
        
        
