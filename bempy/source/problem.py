from .blade_geometry import BladeGeometry
from .results import Result
from .correction_factors import *
from .utils import *
from math import atan, sin, cos, pi
from numpy import isnan, rad2deg
from scipy.optimize import root_scalar
import warnings


class problem:
    def __init__(self, silent_mode=False):
        pass

    def set_parameters(self, rot_speed, wind_speed, rho=1.29):
        if not self.silent_mode:
            print(f"{'#'*10}  PARAMETERS  {'#'*10}")

        self.tip_speed_ratio = self.blade.R * rot_speed / wind_speed
        self.rot_speed = rot_speed
        self.wind_speed = wind_speed
        self.rho = rho

        if not self.silent_mode:
            print(f"rho: {rho}kgm^-3")
            print(f"wind speed: {wind_speed}ms^-1")
            print(f"rotational speed: {rot_speed}ms^-1")
            print(f"tip speed ratio: {self.tip_speed_ratio}")

            print()
            print()

    def add_configuration(self, config_path: str):
        if not self.silent_mode:
            print(f"{'#'*10}  CONFIG  {'#'*10}")

        self.blade.read_configuration(config_path)

        # now initialise the power and thrust element arrays with the correct size
        self.thrust_elements = [0]*(self.blade.number_of_nodes)
        self.torque_elements = [0]*(self.blade.number_of_nodes)
        self.power_elements = [0]*(self.blade.number_of_nodes)

        self.converged_elements = [False]*(self.blade.number_of_nodes)
        self.iterations_elements = [-1]*(self.blade.number_of_nodes)

        if not self.silent_mode:
            print()
            print()

    def apply_ICs(self, axial_IC=-1, tangential_IC=-1):
        if axial_IC == -1:
            axial_IC = [1/3] * self.blade.number_of_nodes

        if tangential_IC == -1:
            tangential_IC = [0] * self.blade.number_of_nodes

        if not self.silent_mode:
            print(f"{'#'*10}  INITIAL CONDITIONS  {'#'*10}")

        if (len(axial_IC) != self.blade.number_of_nodes):
            raise ("Axial ICs wrong size.")
        if (len(tangential_IC) != self.blade.number_of_nodes):
            raise ("Tangential ICs wrong size.")

        self.blade.axial_inductance = axial_IC
        self.blade.tangential_inductance = tangential_IC
        self.blade.calculate_solidity()
        self.blade.calculate_chord_solidity()
        for node in range(self.blade.number_of_nodes):
            self.update_phi(node)
            # print("updating)
        self.blade.update_angle_of_attack()
        self.blade.update_drag_and_lift()
        if not self.silent_mode:
            print("Node:", "Axial IC:", "Tangential IC:")
            for node in range(self.blade.number_of_nodes):
                print(node, round(axial_IC[node], 3),
                      round(tangential_IC[node], 3))
            print()
            print()

    # update the phi angle using the current axial and tangential inducantce factors
    def update_phi(self, specific_node):
        if specific_node in range(self.blade.number_of_nodes):
            axial_inductance = self.blade.axial_inductance[specific_node]
            tangential_inductance = self.blade.tangential_inductance[specific_node]
            radial_distance = self.blade.radial_distances[specific_node]

            phi = self.calculate_phi(
                axial_inductance, tangential_inductance, radial_distance)

            self.blade.phi[specific_node] = phi

        else:
            raise ("ERR: node out of range!")

    def calculate_phi(self, axial_inductance, tangential_inductance, radial_distance):
        temp = self.wind_speed/self.rot_speed
        numer = (1-axial_inductance) * temp
        denom = (1+tangential_inductance) * radial_distance

        phi = atan(numer/denom)

        if isnan(phi):
            raise ("ERR: phi diverged!")

        return phi

    # update the axial and tangential inductance factors using the current phi value

    def update_factors(self, specific_node):
        if specific_node in range(self.blade.number_of_nodes):
            phi = self.blade.phi[specific_node]
            Cl = self.blade.lift_coeff[specific_node]
            Cd = self.blade.drag_coeff[specific_node]
            chord_solidity = self.blade.chord_solidity[specific_node]

            axial_inductance, tangential_inductance = self.calculate_factors(
                phi, Cl, Cd, chord_solidity)

            self.blade.axial_inductance[specific_node] = axial_inductance
            self.blade.tangential_inductance[specific_node] = tangential_inductance
        else:
            raise ("ERR: node out of range!")

    def calculate_factors(self, phi, Cl, Cd, chord_solidity):
        sinphi = sin(phi)
        cosphi = cos(phi)
        Cx = Cl * cosphi + Cd * sinphi
        Cy = Cl * sinphi - Cd * cosphi

        axial_inductance = 1 / (4*sinphi**2 / (chord_solidity * Cx) + 1)
        tangential_inductance = 1 / (4*sinphi*cosphi/(chord_solidity * Cy) - 1)

        return axial_inductance, tangential_inductance

    # calcualtes the power from the tangential and axial inductances
    def calculate_torque(self):
        for i in range(self.blade.number_of_nodes):
            self.torque_elements[i] = self.rho * 4 * pi * self.blade.radial_distances[i]**3 * self.blade.radial_differences[i] * \
                self.blade.tangential_inductance[i] * (
                    1 - self.blade.axial_inductance[i]) * self.rot_speed * self.wind_speed

    def calculate_thrust(self):
        for i in range(self.blade.number_of_nodes):
            self.thrust_elements[i] = self.rho * 4 * pi * self.blade.radial_distances[i] * self.blade.radial_differences[i] * \
                self.blade.axial_inductance[i] * \
                (1 - self.blade.axial_inductance[i]) * self.wind_speed**2

    def get_torque_and_thrust(self):
        return sum(self.torque_elements), sum(self.thrust_elements)

    def apply_tip_loss(self):
        R = self.blade.R
        B = self.blade.B
        tsr = self.tip_speed_ratio
        for i, (r, a, phi) in enumerate(zip(self.blade.radial_distances, self.blade.axial_inductance, self.blade.phi)):
            # factor = prandtl_tip_loss_factor(r, R, B, phi)
            factor = prandtl_tip_loss_factor_approx(r, R, B, a, tsr)
            # factor = xu_sankar_tip_loss(r, R, B, phi)

            self.torque_elements[i] *= factor
            self.thrust_elements[i] *= factor

    def apply_hub_loss(self):
        R_hub = self.blade.R_hub
        B = self.blade.B
        for i, (r, a, phi) in enumerate(zip(self.blade.radial_distances, self.blade.axial_inductance, self.blade.phi)):
            factor = prandtl_hub_loss_factor(r, R_hub, B, phi)

            self.torque_elements[i] *= factor
            self.thrust_elements[i] *= factor

    def calculate_results(self):
        self.calculate_thrust()
        self.calculate_torque()
        (self.torque, self.thrust) = self.get_torque_and_thrust()
        self.power = self.torque * self.rot_speed
        
        
        return self.torque, self.thrust, self.power
    
    def calculate_coefficients(self):
        self.power_coeff = C_power(self.power, self.rho, self.wind_speed, self.blade.R)
        self.thrust_coeff = C_thrust(self.thrust, self.rho, self.wind_speed, self.blade.R)
        self.torque_coeff = C_torque(self.torque, self.rho, self.wind_speed, self.blade.R)
        return self.torque_coeff, self.thrust_coeff, self.power_coeff


class steadyProblem(problem):
    def __init__(self, silent_mode=False):
        self.silent_mode = silent_mode

        self.blade = BladeGeometry(silent_mode)
        self.tip_speed_ratio = 0
        self.rot_speed = 0
        self.wind_speed = 0
        self.rho = 1.29  # density of air (kg m^-3)

        self.err = []

        self.torque_elements = []
        self.thrust_elements = []
        # self.power_elements = [] should add this capability

        self.torque = None
        self.thrust = None
        self.power = None

        self.torque_coeff = None
        self.thrust_coeff = None
        self.power_coeff = None

        self.converged_elements = []
        self.iterations_elements = []

        self.methods = ["iterative", "root-find"]

    # produces a single run until either convergence to the defined tolerance or max interations is met.
    def compute_factors_iterratively(self, tol: float, iter_max: int, show_iterations=False, show_results=True):

        if not self.silent_mode:
            print(f"{'#'*10}  COMPUTE FACTORS (ITERATION)  {'#'*10}")

        for node in range(self.blade.number_of_nodes):
            iter = 0
            converged = False
            while iter < iter_max:
                # change previous factors stored
                self.blade.axial_inductance_prev[node] = self.blade.axial_inductance[node]
                self.blade.tangential_inductance_prev[node] = self.blade.tangential_inductance[node]

                # update the factors
                self.update_factors(node)

                if show_iterations:
                    print("iteration: {iter} err: {err:<10.5f}")

                self.update_phi(node)
                self.blade.update_angle_of_attack(node)
                self.blade.update_drag_and_lift(node)

                err = abs(self.blade.axial_inductance[node] - self.blade.axial_inductance_prev[node])/abs(
                    self.blade.axial_inductance[node])
                if err < tol:
                    converged = True
                    self.converged_elements[node] = True
                    break
                iter += 1

            self.iterations_elements[node] = iter

            if not self.silent_mode:
                print(
                    f"node:{node:<4} converged:{converged:<3} final iteration:{iter:<7} final err:{err:<10.8f}")
        if not self.silent_mode:
            print()
            print()

        if (show_results) and (not self.silent_mode):
            print(f"{'#'*10}  FINAL FACTORS  {'#'*10}")
            for node in range(self.blade.number_of_nodes):
                print(
                    f"{node:<5} a={self.blade.axial_inductance[node]:<10.4f} a'={self.blade.tangential_inductance[node]:<10.4f} phi={rad2deg(self.blade.phi[node]):<7.4f} [deg] twist={rad2deg(self.blade.twist_angles[node]):<7.4f} [deg] alpha={self.blade.angle_of_attack[node]:<7.4f} [deg]")
            print()
            print()

    def get_residual(self, node, phi):
        # update the phi to the new test value
        self.blade.phi[node] = phi

        # update angle of attack and get new drag and lift
        self.blade.update_angle_of_attack(node)
        self.blade.update_drag_and_lift(node)

        radial_distance = self.blade.radial_distances[node]
        Cl = self.blade.lift_coeff[node]
        Cd = self.blade.drag_coeff[node]
        chord_solidity = self.blade.chord_solidity[node]

        axial_inductance, tangential_inductance = self.calculate_factors(
            phi, Cl, Cd, chord_solidity)

        sinphi = sin(phi)
        cosphi = cos(phi)

        temp = self.wind_speed/self.rot_speed

        term1 = sinphi/(1-axial_inductance)
        term2 = cosphi*temp/(radial_distance*(1+tangential_inductance))

        return term1 - term2

    def compute_factors_by_root_finder(self, tol: float, iter_max: int, show_results=True):
        if not self.silent_mode:
            print(f"{'#'*10}  COMPUTE FACTORS (ROOT FINDER)  {'#'*10}")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for node in range(self.blade.number_of_nodes):

                # define lambda function for residual given a node
                def get_residual_node(phi): return self.get_residual(node, phi)
                phi_solution = root_scalar(get_residual_node, method="brentq", bracket=[
                                           0, pi/2], rtol=tol, maxiter=iter_max)

                # now use the root
                self.blade.phi[node] = phi_solution.root
                self.blade.update_angle_of_attack(node)
                self.blade.update_drag_and_lift(node)

                self.update_factors(node)

                self.converged_elements[node] = phi_solution.converged
                self.iterations_elements[node] = phi_solution.iterations

                if not self.silent_mode:
                    print(
                        f"node:{node:<4} converged:{phi_solution.converged:<3} final iteration:{phi_solution.iterations:<7}")

        if not self.silent_mode:
            print()
            print()

        if (show_results) and (not self.silent_mode):
            print(f"{'#'*10}  FINAL FACTORS  {'#'*10}")
            for node in range(self.blade.number_of_nodes):
                print(
                    f"{node:<5} a={self.blade.axial_inductance[node]:<10.4f} a'={self.blade.tangential_inductance[node]:<10.4f} phi={rad2deg(self.blade.phi[node]):<7.4f} [deg] twist={rad2deg(self.blade.twist_angles[node]):<7.4f} [deg] alpha={self.blade.angle_of_attack[node]:<7.4f} [deg]")
            print()
            print()

    def get_methods(self):
        return self.methods()


class UnsteadyProblem(problem):
    def __init__(self, silent_mode=False):
        self.silent_mode = silent_mode

        self.blade = BladeGeometry(silent_mode)

        # constant values for now
        self.tip_speed_ratio = 0
        self.rot_speed = 1.1
        self.wind_speed = 12.1

        self.rho = 1.29  # density of air (kg m^-3)

        self.err = []

        self.torque_elements = []
        self.thrust_elements = []
        # self.power_elements = [] should add this capability

        self.torque = None
        self.thrust = None
        self.power = None

        self.torque_coeff = None
        self.thrust_coeff = None
        self.power_coeff = None

        self.converged_elements = []
        self.iterations_elements = []

        self.methods = ["iterative", "root-find"]

        # additional quantities required for dynamic simulation
        self.t = 0  # initial time
        self.dt = 0.1  # set as 0 for now

        # parameters for the Oye dynamic inflow model
        self.tau_1 = None
        self.tau_2 = None
        self.k = 0.6
  
        self.W_prev = None
        self.W_int_prev = None
        self.W_qs_prev = None
        
        # parameters for Oye dynamic stall model
        self.k_tau = 3
        self.tau_stall = None
        
        self.seperation_func_prev = None
        

        # parameters for original equilibrium simulation
        self.equilibrium_tol = 1E-3
        self.equilibrium_max_iter = 100
        
    # for dynamic stall calculations
        
    

    def update_rot_speed(self, t):
        self.rot_speed = self.rot_speed  # for now this does nothnig to change the dynamics

    def update_wind_speed(self, t, func = None):
        if func == None:
            self.wind_speed = self.wind_speed  # again will do nothing to change the dynamics
        else:
            self.wind_speed = func(t)

    def increase_time_step(self):
        self.t += self.dt

    def filter_induced_velocity(self):
        import numpy as np
        t = self.t
        dt = self.dt
        k = self.k
        wind_speed = self.wind_speed
        rot_speed = self.rot_speed
        R = self.blade.R

        axial_factor = np.array(self.blade.axial_inductance)
        tangential_factor = np.array(self.blade.tangential_inductance)
        radial_distances = np.array(self.blade.radial_distances)

        axial_induced_velocity = axial_factor * wind_speed
        tangential_induced_velocity = tangential_factor*rot_speed*radial_distances

        tau_1 = 1.1*R / ((1-1.3*axial_factor) * wind_speed)
        tau_2 = (0.39 - 0.26*np.power(radial_distances/R, 2)) * tau_1

        tau_1 = np.column_stack((tau_1, tau_1))
        tau_2 = np.column_stack((tau_2, tau_2))

        W_qs_prev = self.W_qs_prev
        W_int_prev = self.W_int_prev
        W_prev = self.W_prev

        # the quasi steady induced flow
        W_qs_cur = np.column_stack(
            (axial_induced_velocity, tangential_induced_velocity))
        H = W_qs_cur + k * tau_1 * (W_qs_cur - W_qs_prev)/dt
        W_int_cur = H + (W_int_prev-H)*np.exp(-dt/tau_1)
        W_cur = W_int_cur + (W_prev - W_int_cur)*np.exp(-dt/tau_2)

        # now update the previous values to the current ones
        self.W_prev = W_cur
        self.W_int_prev = W_int_cur
        self.W_qs_prev = W_qs_cur

    def update_phi_from_previous_induced_vel(self):
        from numpy import arctan, rad2deg
        from math import atan
        wind_speed = self.wind_speed

        for node in range(self.blade.number_of_nodes):
            self.blade.phi[node] = atan((wind_speed - self.W_prev[node][0])
                                        / (self.blade.radial_distances[node]*self.rot_speed + self.W_prev[node][1]))


    def update_inductance_factors_from_previous_induced_vel(self):

        for node in range(self.blade.number_of_nodes):
            self.blade.axial_inductance[node] = self.W_prev[node][0] / \
                self.wind_speed
            self.blade.tangential_inductance[node] = self.W_prev[node][1] / (
                self.rot_speed*self.blade.radial_distances[node])

    def equilibrate_variables(self):
        # here do the intial equilibrium finding that will be used for the first step of the simulation

        # use root-finding method to find equilibrium
        if not self.silent_mode:
            print(f"{'#'*10}  INITIAL EQUILIBRIUM (ROOT FINDER)  {'#'*10}")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for node in range(self.blade.number_of_nodes):

                # define lambda function for residual given a node
                def get_residual_node(phi): return self.get_residual(node, phi)
                phi_solution = root_scalar(get_residual_node, method="brentq", bracket=[
                                           0, pi/2], rtol=self.equilibrium_tol, maxiter=self.equilibrium_max_iter)

                # now use the root
                self.blade.phi[node] = phi_solution.root
                self.blade.update_angle_of_attack(node)
                self.blade.update_drag_and_lift(node)

                self.update_factors(node)

                self.converged_elements[node] = phi_solution.converged
                self.iterations_elements[node] = phi_solution.iterations

                if not self.silent_mode:
                    print(
                        f"node:{node:<4} converged:{phi_solution.converged:<3} final iteration:{phi_solution.iterations:<7}")

        if not self.silent_mode:
            print()
            print()


    def initialise_induced_velocities(self):
        from numpy import array, column_stack
        # now add the initial induction velocities to the arrays
        wind_speed = self.wind_speed
        axial_factor = array(self.blade.axial_inductance)
        tangential_factor = array(self.blade.tangential_inductance)
        radial_distances = array(self.blade.radial_distances)

        axial_induced_velocity = axial_factor * wind_speed
        tangential_induced_velocity = tangential_factor*self.rot_speed*radial_distances

        W = column_stack((axial_induced_velocity, tangential_induced_velocity))
        W = W
        self.W_prev = W
        self.W_int_prev = W
        self.W_qs_prev = W

    def get_residual(self, node, phi):
        # update the phi to the new test value
        self.blade.phi[node] = phi

        # update angle of attack and get new drag and lift
        self.blade.update_angle_of_attack(node)
        self.blade.update_drag_and_lift(node)

        radial_distance = self.blade.radial_distances[node]
        Cl = self.blade.lift_coeff[node]
        Cd = self.blade.drag_coeff[node]
        chord_solidity = self.blade.chord_solidity[node]

        axial_inductance, tangential_inductance = self.calculate_factors(
            phi, Cl, Cd, chord_solidity)

        sinphi = sin(phi)
        cosphi = cos(phi)

        temp = self.wind_speed/self.rot_speed

        term1 = sinphi/(1-axial_inductance)
        term2 = cosphi*temp/(radial_distance*(1+tangential_inductance))

        return term1 - term2
