from BEMsolver import *
import numpy as np
import matplotlib.pyplot as plt



config_path = "/home/dylan/BEMSolver/src/reference_cases/NREL_5MW.case"
problem = Problem()
problem.add_configuration(config_path)

# ICS
axial_IC = [0.33]*17
tangential_IC = [0.1]*17

rated_rpm = 12.1
rated_omega = 2 * pi * rated_rpm / 60
rated_wind_speed = 11.4


problem.set_parameters(rot_speed=rated_omega,wind_speed=rated_wind_speed)
problem.apply_ICs(axial_IC, tangential_IC)
problem.compute_factors(1E-6, 200)
problem.calculate_thrust()
problem.calculate_torque()
problem.apply_tip_loss()
problem.apply_hub_loss()
problem.append_new_results()


print(f"Torque {problem.torque[-1]/1E6} MNm")
print(f"Thrust {problem.thrust[-1]/1E6} MN")
print(f"Power {problem.power[-1]/1E6} MW")


# plt.plot(np.log10(problem.err), label="Error")
# plt.title(f"Residuals: (Tip-Speed-Ratio {problem.tip_speed_ratio:.4f})")
# plt.ylabel("log(Err)")
# plt.xlabel("Iteration Number")
# plt.legend()
# plt.grid()
# plt.show()

