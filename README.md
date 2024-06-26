Warning!: This project is very much still in development meaning the documentation/comments and consistency between modules is poor. This will be addressed soon.

# BEMSolver
bempy, a **B**lade **E**lement **M**omentum theory **PY**thon package for rotor simulations. 
  
```
                                         _                                
                                        | |                               
                                        | |__   ___ _ __ ___  _ __  _   _ 
                                        | '_ \ / _ \ '_ ` _ \| '_ \| | | |
                                        | |_) |  __/ | | | | | |_) | |_| |
                                        |_.__/ \___|_| |_| |_| .__/ \__, |
                                                             | |     __/ |
                                                             |_|    |___/ 
```


# Installation
First clone the repository
```
git clone https://github.com/d-green1958/BEMSolver.git
```
To install using pip run the following command
```
cd BEMSolver
cd bempy
pip install .
```


# Example Usage
```python
# import the package
import bempy

# define a dynamic simulation
sim = bempy.simulation.DynamicSim()

# set the initial conditions
sim.initial_wind_speed = 12.1 #ms^-1
sim.initial_rot_speed = 1.1 #rads^-1

# set the maximum time and time step
sim.t_max = 20
sim.dt = 0.1

# set the turbine data
config_path = bempy.reference_case_path()+"/NREL_5MW.case"
sim.configuration_file = config_path

# tell the simulation to show run time results
sim.show_runtime_results = True

# check the simulation is ready and run!
if sim.check_ready(False):
    sim.run()
```
