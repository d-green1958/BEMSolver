from BEMsolver import *
import numpy as np
import matplotlib.pyplot as plt

def calculate_averages(arr):
    result = []

    for i in range(len(arr)):
        start_index = max(0, i - 2)
        end_index = min(len(arr), i + 3)  # Adjusted to include next 2 values

        subset = arr[start_index:end_index]
        average = sum(subset) / len(subset)
        result.append(average)

    return result



config_path = "/home/dylan/BEMSolver/src/reference_cases/NREL_5MW.case"
problem = Problem()
problem.add_configuration(config_path)

# ICS
axial_IC = [0.3]*17
tangential_IC = [0.8]*17
problem.apply_ICs(axial_IC, tangential_IC, tsr=6)

problem.run(10E-10, 100)
print("AXIAL FACTORS")
print(problem.blade.axial_inductance)

print()
print("TANGENTIAL FACTORS")
print(problem.blade.tangential_inductance)

avg_err = calculate_averages(problem.err)
    

plt.plot(np.log10(problem.err), label="Error")
plt.plot(np.log10(avg_err), label="Average Error")
plt.title(f"Residuals: (Tip-Speed-Ratio {problem.tip_speed_ratio})")
plt.ylabel("log(Err)")
plt.xlabel("Iteration Number")
plt.legend()
plt.grid()
plt.show()
