from math import pi
import matplotlib.pyplot as plt

# power, torque and thrust coefficients
def C_power(power, rho, U, A):
        return 2 * power / (rho * U**3 * A)
    
def C_torque(torque, rho, U, A):
        return 2 * torque / (rho * U**2 * A)
    
def C_thrust(thrust, rho, U, A):
        return 2 * thrust / (rho * U**2 * A)
    
    
# plotting





