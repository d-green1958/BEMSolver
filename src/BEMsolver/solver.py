from .blade_geometry import BladeGeometry
from .correction_factors import *
from math import atan, sin, cos, sqrt, pi
from numpy import isnan, deg2rad, rad2deg


class Problem:
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

        self.torque = []
        self.thrust = []
        self.power = []

        self.torque_coeff = []
        self.thrust_coeff = []
        self.power_coeff = []

    # set the parameters used for solving (rotational speed, wind speed, density of air)

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

        if not self.silent_mode:
            print()
            print()

    def apply_ICs(self, axial_IC, tangential_IC: list[float]):

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
        self.update_phi()
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
    def update_phi(self):
        temp = self.wind_speed/self.rot_speed
        numer = 0
        denom = 0

        for node in range(self.blade.number_of_nodes):
            numer = (1-self.blade.axial_inductance[node])*temp
            denom = (
                1+self.blade.tangential_inductance[node])*self.blade.radial_distances[node]
            self.blade.phi[node] = atan(numer/denom)

            if isnan(self.blade.phi[node]):
                print("ERR: DIVERGENCE")
                raise ("ERR: Iteration has diverged!")

    # update the axial and tangential inductance factors using the current phi value
    def update_factors(self):
        B = self.blade.B
        Cx = 0
        Cy = 0

        sinphi = 0
        cosphi = 0

        for node, (phi, Cl, Cd, sigma) in enumerate(zip(self.blade.phi,
                                                        self.blade.lift_coeff,
                                                        self.blade.drag_coeff,
                                                        self.blade.chord_solidity)):
            sinphi = sin(phi)
            cosphi = cos(phi)
            Cx = Cl * cosphi + Cd * sinphi
            Cy = Cl * sinphi - Cd * cosphi

            self.blade.axial_inductance_prev[node] = self.blade.axial_inductance[node]
            self.blade.tangential_inductance_prev[node] = self.blade.tangential_inductance[node]

            self.blade.axial_inductance[node] = 1 / \
                (4*sinphi**2 / (sigma * Cx) + 1)
            self.blade.tangential_inductance[node] = 1 / \
                (4*sinphi*cosphi/(sigma * Cy) - 1)

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
        

    def append_new_results(self):
        # thrust = self.calculate_thrust()
        # torque = self.calculate_torque()

        (torque, thrust) = self.get_torque_and_thrust()

        power = torque * self.rot_speed

        self.thrust.append(thrust)
        self.torque.append(torque)
        self.power.append(power)

    # returns L_2 error squared

    def find_err(self):
        sum = 0
        for node in range(self.blade.number_of_nodes):
            sum += (self.blade.axial_inductance[node] - self.blade.axial_inductance_prev[node])**2 + (
                self.blade.tangential_inductance[node] - self.blade.tangential_inductance_prev[node])**2
        return sum

    # produces a single run until either convergence to the defined tolerance or max interations is met.
    def compute_factors(self, tol: float, iter_max: int, show_iterations=False, show_results=True):
        if not self.silent_mode:
            print(f"{'#'*10}  COMPUTE FACTORS  {'#'*10}")

        converged = False
        iter = 0
        while iter < iter_max:
            self.update_factors()
            self.err.append(sqrt(self.find_err()))
            if show_iterations:
                print("iteration: ", iter, "err: ", self.err[-1])

            self.update_phi()
            self.blade.update_angle_of_attack()
            self.blade.update_drag_and_lift()

            if (isnan(self.err[-1])):
                break

            if (self.err[-1] < tol):
                break
            iter += 1

        if not self.silent_mode:
            print(f"converged:{converged} final iteration:{iter}")
            print()
            print()

        if (show_results) and (not self.silent_mode):
            print(f"Final Values:")
            for node in range(self.blade.number_of_nodes):
                print(
                    f"{node:<5} a={self.blade.axial_inductance[node]:<10.4f} a'={self.blade.tangential_inductance[node]:<10.4f} phi={rad2deg(self.blade.phi[node]):<10.4f} [deg]")
            print()
            print()

        # produced multiple runs over provided parameters (tip speeds).
        def multiple_run(self, tol: float, iter_max: int, rot_speeds: list[float], wind_speeds: list[float], axial_IC, tangnetial_IC, show_iterations=True):
            print("FINISH THIS")
            if not self.silent_mode():
                print(f"{'#'*10}  MULTIPLE SOLVE  {'#'*10}")

            # change parameters
            # apply the ICs
            # run until convergence
            # calcualte power and C_p
            # save to
