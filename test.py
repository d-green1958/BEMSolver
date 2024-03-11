from math import pi

import PyBEM as BEM


# ICS
axial_IC = [0.33]*17
tangential_IC = [0.1]*17

# parameters
rated_rpm = 12.1
rated_omega = 2 * pi * rated_rpm / 60
rated_wind_speed = 11.4

# config file
config_path = "/home/dylan/BEMSolver/PyBEM/reference_cases/NREL_5MW.case"

# run
sol = BEM.single_run(config_path, rated_wind_speed, rated_omega, use_hub_loss=True, use_tip_loss=True, tol=1E-6, max_iterations=1000)
sol1 = BEM.single_run(config_path, rated_wind_speed, rated_omega, use_hub_loss=True, use_tip_loss=True, tol=1E-6, max_iterations=1000)

multisol = BEM.parametric_run(config_path, rated_wind_speed-10, rated_wind_speed+2, 5,
                     -1+rated_omega, rated_omega+1, 5,
                     use_tip_loss=True, use_hub_loss=True,
                     silent_mode=True, show_runtime_results=True, 
                     method="root-find", tol=1E-3)


# # Plot
# BEM.plotBladeLoads(sol1)
# BEM.plotBladeLoads(sol)

# BEM.plotCPowerTsr(multisol)
# BEM.plotCTorqueTsr(multisol)
# BEM.plotCThrustTsr(multisol)


# multisol1 = parametric_run(config_path, rated_wind_speed, rated_wind_speed, 1,
#                      -1+rated_omega, rated_omega+1, 50,
#                      use_tip_loss=False, use_hub_loss=False,
#                      silent_mode=True, show_runtime_results=True, 
#                      method="root-find", tol=1E-3)

# plotCPowerTsr(multisol, add_to_last_plot=False)
# plotCPowerTsr(multisol1, add_to_last_plot=True)

# plotCTorqueTsr(multisol, add_to_last_plot=False)
# plotCTorqueTsr(multisol1, add_to_last_plot=True)

# plotCThrustTsr(multisol, add_to_last_plot=False)
# plotCThrustTsr(multisol1, add_to_last_plot=True)


# BEM.show_plots()


   
