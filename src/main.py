from BEMsolver import *

# ICS
axial_IC = [0.33]*17
tangential_IC = [0.1]*17

rated_rpm = 12.1
rated_omega = 2 * pi * rated_rpm / 60
rated_wind_speed = 11.4

config_path = "/home/dylan/BEMSolver/src/reference_cases/NREL_5MW.case"

single_run(config_path, rated_wind_speed, rated_omega, use_hub_loss=True, use_tip_loss=True)
# parametric_run(config_path, 10, 11, 2, 1, 1.1, 2)