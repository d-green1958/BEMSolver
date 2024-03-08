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


sol = parametric_run(config_path, 11.9, 12.4, 20, -1+rated_omega, rated_omega+2, 50,
                     use_tip_loss=False, use_hub_loss=False, silent_mode=True, show_runtime_results=True, 
                     method="root-find", tol=1E-3)


plotCPowerTsr(sol)
plotCTorqueTsr(sol)
plotCThrustTsr(sol)
plt.show()


# max_val = 0
# for i in sol:
#     for j in i:
#         # print(j.power_coefficient)
#         # print(j.sectional_convergence)
#         if j.power_coefficient > max_val:
#             max_val = j.power_coefficient
   
# # TSR = []
# CP = []         
# for i in sol:
#     for j in i:
#         TSR.append(j.tip_speed_ratio)
#         CP.append(j.power_coefficient)
                     

# print(max_val)


    

# plt.plot(TSR, CP, marker="+")
# plt.grid()
# plt.hlines([16/27],linestyles="--", xmin=0, xmax=max(TSR),colors="red", label="Betz")
# plt.legend()
# plt.show()
    

