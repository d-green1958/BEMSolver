from .blade_geometry import BladeGeometry
from math import atan,sin,cos,sqrt, pi
from numpy import isnan
    
class Problem:
    def __init__(self, silent_mode=False):
        self.silent_mode = silent_mode
        
        self.blade = BladeGeometry(silent_mode)
        self.tip_speed_ratio = 0
        self.rot_speed = 0
        self.wind_speed = 0
        self.rho = 1.29 # density of air (kg m^-3)

        self.err = []
        
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
        
        if not self.silent_mode:
            print()
            print()
        
    def apply_ICs(self, axial_IC, tangential_IC: list[float]):
        if not self.silent_mode:
            print(f"{'#'*10}  INITIAL CONDITIONS  {'#'*10}")
        
        if (len(axial_IC) != self.blade.number_of_nodes):
            raise("Axial ICs wrong size.")
        if (len(tangential_IC) != self.blade.number_of_nodes):
            raise("Tangential ICs wrong size.")
        
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
                print(node, round(axial_IC[node],3), round(tangential_IC[node],3))
            print()
            print()
    
    # update the phi angle using the current axial and tangential inducantce factors
    def update_phi(self):
        temp = self.blade.R/self.tip_speed_ratio
        numer = 0
        denom = 0

        for node in range(self.blade.number_of_nodes):
            numer = (1-self.blade.axial_inductance[node])*temp
            denom = (1+self.blade.tangential_inductance[node])*self.blade.radial_distances[node]
            self.blade.phi[node] = atan(numer/denom)
            
            if isnan(self.blade.phi[node]):
                print("ERR: DIVERGENCE")
                # print("NODE", node)
                # print("NUMERATOR", (1-self.blade.axial_inductance[node]))
                # print("TEMP", temp)
                # print("DENOMINATOR", (1+self.blade.tangential_inductance[node]))
                # print("DR", self.blade.radial_distances[node])

    # update the axial and tangential inductance factors using the current phi value
    def update_factors(self):
        B = self.blade.B
        Cx = 0
        Cy = 0
        
        sinphi = 0
        cosphi = 0
                
        for node,(phi, Cl, Cd, sigma) in enumerate(zip(self.blade.phi,
                                               self.blade.lift_coeff,
                                               self.blade.drag_coeff,
                                               self.blade.chord_solidity)):
            sinphi = sin(phi)
            cosphi = cos(phi)
            Cx = Cl * cosphi + Cd * sinphi
            Cy = Cl * sinphi - Cd * cosphi
            
            self.blade.axial_inductance_prev[node] = self.blade.axial_inductance[node]
            self.blade.tangential_inductance_prev[node] = self.blade.tangential_inductance[node]

            self.blade.axial_inductance[node] = 1/(4*sinphi**2 / (sigma * Cx) + 1)
            self.blade.tangential_inductance[node] = 1/(4*sinphi*cosphi/(sigma * Cy) - 1)
         
    # calcualtes the power from the tangential and axial inductances
    def calculate_torque(self, U_inf, rot_speed):
        temp = 0
        for i in range(self.blade.number_of_nodes):
            temp += self.rho* 4 * pi * self.blade.radial_distances[i]**3 * self.blade.radial_differences[i] * self.blade.tangential_inductance[i] * (1 - self.blade.axial_inductance[i]) * rot_speed
        return temp
    
    # returns L_2 error squared
    def find_err(self):
        sum = 0
        for node in range(self.blade.number_of_nodes):
            sum += (self.blade.axial_inductance[node] - self.blade.axial_inductance_prev[node])**2 +(self.blade.tangential_inductance[node] - self.blade.tangential_inductance_prev[node])**2
        return sum
            
    # produces a single run until either convergence to the defined tolerance or max interations is met.
    def single_run(self, tol: float, iter_max: int):    
        if not self.silent_mode:
            print(f"{'#'*10}  SINGLE SOLVE  {'#'*10}")
            
        converged = False
        iter = 0
        print()
        while iter<iter_max:
            self.update_factors()
            self.err.append(sqrt(self.find_err()))
            if not self.silent_mode:
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
            
        
        # produced multiple runs over provided parameters (tip speeds).
        def multiple_run(self, tol:float, iter_max: int, rot_speeds: list[float], wind_speeds: list[float], axial_IC, tangnetial_IC):
            print("FINISH THIS")
            if not self.silent_mode():
                print(f"{'#'*10}  MULTIPLE SOLVE  {'#'*10}")
            
            # change parameters
            # apply the ICs
            # run until convergence
            # calcualte power and C_p
            # save to 
                
            
            
            
        
    