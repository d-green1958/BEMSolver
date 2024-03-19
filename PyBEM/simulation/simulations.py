# from source import steadyProblem, Result, C_power, C_thrust, C_torque
# from numpy import ndindex, meshgrid, linspace, ones


class steadySim():
    
    def __init__(self):
        return
            
    def single_run(self, configuration_file, wind_speed, rot_speed, axial_initial = -1, tangential_initial = -1,
                   tol = 1E-5, max_iterations = 100, silent_mode = True, use_tip_loss = False,
                   use_hub_loss = False, method = "root-find"):
        
        from source import steadyProblem, Result, C_power, C_thrust, C_torque
        from numpy import ndindex, meshgrid, linspace, ones
        
        problem = steadyProblem(silent_mode)
        
        methods = problem.methods
        
        if method not in methods:
            raise ValueError("ERR: invalid method!")
        
        problem.add_configuration(config_path=configuration_file)
        
        if axial_initial == -1:
            axial_initial = [1/3]*problem.blade.number_of_nodes
        if tangential_initial == -1:
            tangential_initial = [0]*problem.blade.number_of_nodes
        
        problem.set_parameters(rot_speed=rot_speed,wind_speed=wind_speed)
        problem.apply_ICs(axial_IC=axial_initial, tangential_IC=tangential_initial)
        if method == "iterative":
            problem.compute_factors_iterratively(tol=tol, iter_max=max_iterations)
        if method == "root-find":
            problem.compute_factors_by_root_finder(tol=tol, iter_max=max_iterations)

        # now compute results
        problem.calculate_thrust()
        problem.calculate_torque()
        
        # apply correction factors
        if use_tip_loss:
            problem.apply_tip_loss()
        if use_hub_loss:
            problem.apply_hub_loss()
    
        # compute final values
        torque, thrust, power = problem.calculate_results()    
        
        if not silent_mode:
            print(f"{'#'*10}  FINAL RESULTS  {'#'*10}")
            print(f"Torque {torque/1E6:<10.5f} MNm")
            print(f"Thrust {thrust/1E6:<10.5f} MN")
            print(f"Power  {power/1E6:<10.5f} MW")
        
        # create a results object to return
        results = Result() 
        results.wind_speed = wind_speed
        results.rot_speed = rot_speed
        results.tip_speed_ratio = problem.tip_speed_ratio
        
        results.sectional_axial_inductance_factor = problem.blade.axial_inductance
        results.sectional_tangential_inductance_factor = problem.blade.tangential_inductance
        results.sectional_phi = problem.blade.phi
        
        results.sectional_torque = problem.torque_elements
        results.sectional_thrust = problem.thrust_elements
        results.sectional_positions = problem.blade.radial_distances
        results.sectional_chord_lengths = problem.blade.chord_lengths
        
        results.power = power
        results.thrust = thrust
        results.torque = torque
        
        results.sectional_iterations = problem.iterations_elements
        results.sectional_convergence = problem.converged_elements
        
        
        return results



        
    def parametric_run(configuration_file, wind_speed_start, wind_speed_end, wind_speed_nodes,
                    rot_speed_start, rot_speed_end, rot_speed_nodes,
                    axial_initial=-1, tangential_initial=-1,
                    tol = 1E-5, max_iterations = 100, silent_mode = True, use_tip_loss = True,
                    use_hub_loss = True, method="root-find", show_runtime_results = False,
                    rho = 1.29, B=3): 
        
        from source import steadyProblem, Result, C_power, C_thrust, C_torque
        from numpy import ndindex, meshgrid, linspace, ones
        
        problem = steadyProblem(silent_mode=silent_mode)
        methods = problem.methods
        
        if method not in methods:
            raise ValueError("ERR: invalid method!")   

        problem.add_configuration(config_path=configuration_file)
        problem.blade.B = B
        problem.rho = rho
        
        if axial_initial == -1:
            axial_initial = [1/3]*problem.blade.number_of_nodes
        if tangential_initial == -1:
            tangential_initial = [0]*problem.blade.number_of_nodes
            
        wind_speeds = linspace(wind_speed_start, wind_speed_end, wind_speed_nodes)
        rot_speeds = linspace(rot_speed_start, rot_speed_end, rot_speed_nodes)
        
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
            
            C_T = C_thrust(thrust, rho, wind_speed, problem.blade.R)
            C_P = C_power(power, rho, wind_speed, problem.blade.R)
            C_Q = C_torque(torque, rho, wind_speed, problem.blade.R)
            
            results.torque_coefficient = C_Q
            results.power_coefficient = C_P
            results.thrust_coefficient = C_T
            
            if show_runtime_results:
                print(f"Power:{power/1E6:<8.6}MW Torque:{torque/1E6:<8.6}MNm Thrust:{thrust/1E6:<8.6}MN")
                print(f"C_P:{C_P:<8.6} C_Q:{C_Q:<8.6} C_T:{C_T:<8.6}")
                print()
            
            results_array[wind_speed_index][rot_speed_index] = results
            
        return results_array
    
    
    
    

            
