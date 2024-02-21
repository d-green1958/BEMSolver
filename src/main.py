from BEMsolver import *
import matplotlib.pyplot as plt

# ICS
axial_IC = [0.33]*17
tangential_IC = [0.1]*17

rated_rpm = 12.1
rated_omega = 2 * pi * rated_rpm / 60
rated_wind_speed = 11.4

config_path = "/home/dylan/BEMSolver/src/reference_cases/NREL_5MW.case"

single_run(config_path, rated_wind_speed, rated_omega, use_hub_loss=True, use_tip_loss=True, tol=1E-6, max_iterations=1000)
# single_run(config_path, 13, 1.4, use_hub_loss=True, use_tip_loss=True, tol=1E-6, max_iterations=1000)

# torque, thrust, power, wind_speed, rot_speed = parametric_run(config_path, 9, 13, 20, 0.9, 1.4, 10)




