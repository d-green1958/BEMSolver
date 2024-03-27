
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
        
        
    def __str__(self) -> str:
        return (f"Results:\n"
                f"Wind Speed {self.wind_speed}\n"
                f"Rotational Speed {self.rot_speed}\n"
                f"Power {self.power/1E6:<10.5f} MW ---> C_P {self.power_coefficient:<10.5f}\n"
                f"Torque {self.torque/1E6:<10.5f} MW ---> C_Q {self.torque_coefficient:<10.5f}\n"
                f"Thrust {self.thrust/1E6:<10.5f} MW ---> C_T {self.thrust_coefficient:<10.5f}\n")


class DynamicResult(Result):
    def __init__(self) -> None:
        
        # simulation parameters
        self.wind_speed = []
        self.rot_speed = []
        self.tip_speed_ratio = []
        self.rotor_radius = None
        
        # blade data
        self.sectional_positions = []
        self.sectional_chord_lengths = []
        self.sectional_lengths = []
        
        # sectional results
        self.sectional_torque = []
        self.sectional_thrust = []
        
        # integrated results
        self.power = []
        self.thrust = []
        self.torque = []
        
        # inductance factors and inflow angle
        self.sectional_axial_inductance_factor = []
        self.sectional_tangential_inductance_factor = []
        self.sectional_phi = []
        self.sectional_angle_of_attack = []
        
        # simulation performance
        self.sectional_convergence = []
        self.sectional_iterations = []
        self.sectional_final_error = []
        
        # performance coefficients
        self.power_coefficient = []
        self.torque_coefficient = []
        self.thrust_coefficient = []
        
        # lift and drag coefficients
        self.sectional_lift_coefficient = []
        self.sectional_drag_coefficient = []
        
        self.time = []
        
    def add_data(self, problem):
        """
        function used to quickly add data from a dynamic simulation to the existing results

        Args:
            problem (UnsteadyProblem): the UnsteadyProblem class containing the results
        """
        self.time.append(problem.t)
        
        self.wind_speed.append(problem.wind_speed)
        self.rot_speed.append(problem.rot_speed)
        
        self.sectional_tangential_inductance_factor.append(problem.blade.tangential_inductance.copy())
        self.sectional_axial_inductance_factor.append(problem.blade.axial_inductance.copy())
        self.sectional_phi.append(problem.blade.phi.copy())
        
        self.sectional_thrust.append(problem.thrust_elements.copy())
        self.sectional_torque.append(problem.torque_elements.copy())
        
        self.power.append(problem.power)
        self.power_coefficient.append(problem.power_coeff)
        self.thrust.append(problem.thrust)
        self.thrust_coefficient.append(problem.thrust_coeff)
        self.torque.append(problem.torque)
        self.torque_coefficient.append(problem.torque_coeff)
        
        self.sectional_angle_of_attack.append(problem.blade.angle_of_attack.copy())
        self.sectional_lift_coefficient.append(problem.blade.lift_coeff.copy())
        self.sectional_drag_coefficient.append(problem.blade.drag_coeff.copy())