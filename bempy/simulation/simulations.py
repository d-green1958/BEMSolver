# from source import steadyProblem, Result, C_power, C_thrust, C_torque
# from numpy import ndindex, meshgrid, linspace, ones


class SteadySim():
    
    def __init__(self):
        # inflow velocity parameter range
        self.wind_speeds = None
        self.rot_speeds = None
        
        # default parameters
        self.rho = 1.29
        self.B = 3
        
        # simulation parameters
        self.tol = 1E-5
        self.max_iterations = 100
        self.method = "root-find"
        
        # initial conditions for simulation
        self.axial_initial=-1
        self.tangential_initial=-1
        
        # outstream
        self.show_runtime_results = False
        self.silent_mode = True
        
        # correction factors
        self.use_tip_loss = True
        self.use_hub_loss = True
        
        # inputs
        self.configuration_file = None
        
        # outputs
        self.results = None
        return
    
    
    def show_parameters(self):
        pass # here add a summary of the simulation we are going to run.    
    
    def run(self): 
        from source import steadyProblem, Result, C_power, C_thrust, C_torque
        from numpy import ndindex, meshgrid, linspace, ones
        
        # not efficient but will work for now (oops!)
        tol = self.tol
        max_iterations = self.max_iterations
        axial_initial = self.axial_initial
        tangential_initial = self.tangential_initial
        show_runtime_results = self.show_runtime_results
        method = self.method
        silent_mode = self.silent_mode
        use_tip_loss = self.use_tip_loss
        use_hub_loss = self.use_hub_loss
        configuration_file = self.configuration_file
        wind_speeds = self.wind_speeds
        rot_speeds = self.rot_speeds
        
        problem = steadyProblem(silent_mode=silent_mode)
        methods = problem.methods
        
        if method not in methods:
            raise ValueError("ERR: invalid method!")   

        problem.add_configuration(config_path=configuration_file)
        problem.blade.B = self.B
        problem.rho = self.rho
            
        # wind_speeds = linspace(wind_speed_start, wind_speed_end, wind_speed_nodes)
        # rot_speeds = linspace(rot_speed_start, rot_speed_end, rot_speed_nodes)
        wind_speed_nodes = len(wind_speeds)
        rot_speed_nodes = len(rot_speeds)
        
        wind_speeds, rot_speeds = meshgrid(wind_speeds, rot_speeds)
        
        powers = ones(wind_speeds.shape)
        thrusts = powers
        torques = powers
        
        obj = Result()
        # results_array = [[obj]*wind_speed_nodes]*rot_speed_nodes
        results_array = []
        
        
        for wind_speed_index in range(len(wind_speeds)):
            row = []
            for rot_speed_index in range(len(rot_speeds[wind_speed_index])):
                obj = Result()
                row.append(obj)
            results_array.append(row)
        
        
        counter = 0
        total_counter = wind_speed_nodes*rot_speed_nodes
        for wind_speed_index, rot_speed_index in ndindex(wind_speeds.shape):
        
            counter += 1
            
            wind_speed = wind_speeds[wind_speed_index][rot_speed_index]
            rot_speed = rot_speeds[wind_speed_index][rot_speed_index]
            
            # print(wind_speed)
            # print(rot_speed)
            
            problem.wind_speed = wind_speed
            problem.rot_speed = rot_speed    
                
            tsr = problem.blade.R * rot_speed / wind_speed

            if show_runtime_results:
                print(f"ITERATION NUM:{counter}")
                print(f"PROGRESS:{100*counter/total_counter:.5}%")
                print(f"parameters: wind speed:{wind_speed:<6.5}, rot_speed:{rot_speed:<6.5}, tsr:{tsr:<6.5}")
            
            
            problem.set_parameters(wind_speed=wind_speed, rot_speed=rot_speed)
            problem.apply_ICs(axial_IC=axial_initial, tangential_IC=tangential_initial)
            
            if method == "iterative":
                problem.compute_factors_iterratively(tol=tol, iter_max=max_iterations)
            if method == "root-find":
                problem.compute_factors_by_root_finder(tol=tol, iter_max=max_iterations)


            problem.calculate_thrust()
            problem.calculate_torque()
            
            # apply correction factors
            if use_tip_loss:
                problem.apply_tip_loss()
            if use_hub_loss:
                problem.apply_hub_loss()
            
            torque, thrust, power = problem.calculate_results()
            
            
            results = Result() 
            results.wind_speed = wind_speed
            results.rot_speed = rot_speed
            results.tip_speed_ratio = problem.tip_speed_ratio
            
            results.sectional_axial_inductance_factor = problem.blade.axial_inductance
            results.sectional_tangential_inductance_factor = problem.blade.tangential_inductance
            results.sectional_phi = problem.blade.phi
            
            results.sectional_torque = problem.torque_elements
            results.sectional_thrust = problem.thrust_elements
            
            results.power = power
            results.thrust = thrust
            results.torque = torque
            
            results.rotor_radius = problem.blade.R
            
            results.sectional_iterations = problem.iterations_elements
            results.sectional_convergence = problem.converged_elements
            
            C_T = C_thrust(thrust, self.rho, wind_speed, problem.blade.R)
            C_P = C_power(power, self.rho, wind_speed, problem.blade.R)
            C_Q = C_torque(torque, self.rho, wind_speed, problem.blade.R)
            
            results.torque_coefficient = C_Q
            results.power_coefficient = C_P
            results.thrust_coefficient = C_T
            
            if show_runtime_results:
                print(f"Power:{power/1E6:<8.6}MW Torque:{torque/1E6:<8.6}MNm Thrust:{thrust/1E6:<8.6}MN")
                print(f"C_P:{C_P:<8.6} C_Q:{C_Q:<8.6} C_T:{C_T:<8.6}")
                print()
                
            elif not silent_mode:
                print(f"{'#'*10}  FINAL RESULTS  {'#'*10}")
                print(f"Torque {torque/1E6:<10.5f} MNm ---> C_Q {C_Q:<10.5f}")
                print(f"Thrust {thrust/1E6:<10.5f} MN  ---> C_T {C_T:<10.5f}")
                print(f"Power  {power/1E6:<10.5f} MW  ---> C_P {C_P:<10.5f}")
                print()
                print()
        
            
            
            results_array[wind_speed_index][rot_speed_index] = results
            self.results = results_array
            
        return results_array
    
    
    
    

            
class DynamicSim():
    def __init__(self):
        # inflow velocity parameter range
        # self.wind_speeds = None  # these need to be able to be a functino of time 
        # self.rot_speeds = None

        # default parameters
        self.rho = 1.29
        self.B = 3

        # simulation parameters
        self.tol = 1E-5
        self.max_iterations = 100
        self.method = "root-find"
        self.t_max = 3
        self.dt = 0.2

        # initial conditions for simulation
        self.axial_initial=-1
        self.tangential_initial=-1

        # outstream
        self.show_runtime_results = False
        self.silent_mode = True

        # correction factors
        self.use_tip_loss = True
        self.use_hub_loss = True

        # inputs
        self.configuration_file = None

        # outputs
        self.results = None
        return
    
    
    def run(self):
        from source import UnsteadyProblem
        # 0: configuration
        # 0.1: intial conditions
        # 0.2: find equilibrium for initla conditions
        # 0.3: calculate the induction velocity
        
        # 1: next time step
        # 2: compute next a,a',phi
        # 3: relax quasisteady induction velocity to find induction velocity
        # 4: used induced velocities to find phi, a and a'
        
        # not efficient but will work for now (oops!)
        tol = self.tol
        max_iterations = self.max_iterations
        axial_initial = self.axial_initial
        tangential_initial = self.tangential_initial
        show_runtime_results = self.show_runtime_results
        method = self.method
        silent_mode = self.silent_mode
        use_tip_loss = self.use_tip_loss
        use_hub_loss = self.use_hub_loss
        configuration_file = self.configuration_file
        t_max = self.t_max
        dt = self.dt
        wind_speeds = self.wind_speeds
        rot_speeds = self.rot_speeds
        
        problem = UnsteadyProblem(silent_mode)
        
        problem.set_parameters(wind_speed=wind_speeds[0], rot_speed=rot_speeds[0]) # THIS WILL NEED TO BE CHANGED
        
        # make sure problem has right parameters
        problem.dt = dt
        problem.blade.B = self.B
        problem.rho = self.rho
        
        # add the configuration file
        problem.add_configuration(configuration_file)
        problem.apply_ICs(axial_initial,tangential_initial)

        # find equilibrium values for t=0
        problem.equilibrate_variables()
        problem.initialise_induced_velocities()
        
        ax_for_plot = []
        ax_for_plot.append(problem.blade.axial_inductance[10])
        
        while problem.t < t_max:
            print("TIME", problem.t)
            
            # now increase time step and adjust wind speed and rot speed
            problem.increase_time_step()
            problem.update_wind_speed(problem.t)
            problem.update_rot_speed(problem.t)
            
            # now calculate the new a,a' and phi from previous W
            # problem.update_phi_from_previous_induced_vel()
            # for node in range(problem.blade.number_of_nodes):
            #     problem.update_phi(node)
            #     problem.blade.update_angle_of_attack(node)
            #     problem.blade.update_drag_and_lift(node)
            #     problem.update_factors(node)
            
            
            problem.equilibrate_variables()
            
            # use new a,a' and phi to calculate new W and filter it
            problem.filter_induced_velocity()
            
            
            problem.update_phi_from_previous_induced_vel()
            problem.update_inductance_factors_from_previous_induced_vel()
        
            ax_for_plot.append(problem.blade.axial_inductance[8])
        
        return ax_for_plot
