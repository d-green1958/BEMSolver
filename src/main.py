from BEMsolver import *
import matplotlib.pyplot as plt
from numpy import amin

# ICS
axial_IC = [0.33]*17
tangential_IC = [0.1]*17

rated_rpm = 12.1
rated_omega = 2 * pi * rated_rpm / 60
rated_wind_speed = 11.4

config_path = "/home/dylan/BEMSolver/src/reference_cases/NREL_5MW.case"

# sol = single_run(config_path, rated_wind_speed, rated_omega, use_hub_loss=True, use_tip_loss=True, tol=1E-6, max_iterations=1000)
# sol = single_run(config_path, 13, 1.4, use_hub_loss=True, use_tip_loss=True, tol=1E-6, max_iterations=1000)


sol = parametric_run(config_path, 10, 12, 100, -0.2+rated_omega, 0.2+rated_omega, 100,
                     use_tip_loss=True, use_hub_loss=True, silent_mode=True, show_runtime_results=True)

max_val = 0
for i in sol:
    for j in i:
        # print(j.power_coefficient)
        if j.power_coefficient > max_val:
            max_val = j.power_coefficient

print(max_val)




