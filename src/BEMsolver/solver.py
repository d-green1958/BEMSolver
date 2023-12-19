from .blade_geometry import BladeGeometry
from math import atan,sin,cos,sqrt
from numpy import arctan, isnan
    
class Problem:
    def __init__(self):
        self.blade = BladeGeometry()
        self.tip_speed_ratio = 0
        
        self.err = []
    
    def add_configuration(self, config_path: str):
        boarder = "#"*10
        print(boarder + "  CONFIG  " + boarder)
        self.blade.read_configuration(config_path)
        print()
        print()
        
      
    def apply_ICs(self, axial_IC, tangential_IC: list[float], tsr):
        boarder = "#"*10
        print(boarder + "  INITIAL CONDITIONS  " + boarder)
        self.tip_speed_ratio = tsr
        print("Tip speed ratio:  " + str(tsr))
        
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
        
        print("Node:", "Axial IC:", "Tangential IC:")
        for node in range(self.blade.number_of_nodes):
            print(node, round(axial_IC[node],3), round(tangential_IC[node],3))
        
        print()
        print()
        
    def update_phi(self):
        temp = self.blade.R/self.tip_speed_ratio
        numer = 0
        denom = 0

        for node in range(self.blade.number_of_nodes):
            numer = (1-self.blade.axial_inductance[node])*temp
            denom = (1+self.blade.tangential_inductance[node])*self.blade.radial_distances[node]
            # self.blade.phi[node] = atan((1-self.blade.axial_inductance[node])*temp/((1+self.blade.tangential_inductance[node])*self.blade.radial_distances[node]))
            self.blade.phi[node] = atan(numer/denom)
            # self.blade.phi[node] = atan(self.blade.tangential_inductance[node]*self.blade.radial_distances[node]*self.tip_speed_ratio/(self.blade.axial_inductance[node]*self.blade.R))
            
            
            if isnan(self.blade.phi[node]):
                print("NODE", node)
                print("NUMERATOR", (1-self.blade.axial_inductance[node]))
                print("TEMP", temp)
                print("DENOMINATOR", (1+self.blade.tangential_inductance[node]))
                print("DR", self.blade.radial_distances[node])

        
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
            
            
    
    # returns L_2 error squared
    def find_err(self):
        sum = 0
        for node in range(self.blade.number_of_nodes):
            sum += (self.blade.axial_inductance[node] - self.blade.axial_inductance_prev[node])**2 +(self.blade.tangential_inductance[node] - self.blade.tangential_inductance_prev[node])**2
        return sum
            
         
    
    def run(self, tol: float, iter_max: int):
        boarder = "#"*10
        print(boarder + "  SOLVE  " + boarder)
  
        iter = 0
        print()
        while iter<iter_max:
            self.update_factors()
            self.err.append(sqrt(self.find_err()))
            print("iteration: ", iter, "err: ", self.err[-1])
            
            self.update_phi()
            self.blade.update_angle_of_attack()
            self.blade.update_drag_and_lift()
            
            if (isnan(self.err[-1])):
                break

            if (self.err[-1] < tol):
                print("converged.")
                break
            
            
            iter += 1
        
        print()
        print()
    