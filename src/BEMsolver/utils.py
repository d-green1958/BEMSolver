from math import pi, sqrt
import matplotlib.pyplot as plt

# power, torque and thrust coefficients
def C_power(power, rho, U, R):
        return 2 * power / (rho * U**3 * pi * R**2)
    
def C_torque(torque, rho, U, R):
        return 2 * torque / (rho * U**2 * pi * R**3)
    
def C_thrust(thrust, rho, U, R):
        return 2 * thrust / (rho * U**2 * pi * R**2)
    
    
# plotting





