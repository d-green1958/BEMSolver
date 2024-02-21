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

# single_run(config_path, rated_wind_speed, rated_omega, use_hub_loss=True, use_tip_loss=True, tol=1E-6, max_iterations=1000)
# single_run(config_path, 13, 1.4, use_hub_loss=True, use_tip_loss=True, tol=1E-6, max_iterations=1000)

torque, thrust, power, wind_speeds, rot_speeds = parametric_run(config_path, 4, 25, 100, rated_omega-0.3, rated_omega+0.3, 5, use_tip_loss=True, use_hub_loss=True)

CP = ones(power.shape)
tsr = 63.5*rot_speeds/wind_speeds

for wind_speed_index, rot_speed_index in ndindex(wind_speeds.shape):
    CP[wind_speed_index][rot_speed_index] = C_power(power[wind_speed_index][rot_speed_index], 1.025, wind_speeds[wind_speed_index][rot_speed_index], 65)

for row in range(len(wind_speeds)):
    plt.plot(tsr[row], CP[row], label=rot_speeds[row][0])

plt.ylim([0, 1])
plt.hlines(16/27, xmin=0, xmax=25, linestyles="--", colors="black", label="Betz")
plt.xlim([amin(tsr)-2,amax(tsr)+2])
plt.grid()
plt.xlabel("TSR")
plt.ylabel("C_P")
plt.legend(title="rot speed")
plt.show()




