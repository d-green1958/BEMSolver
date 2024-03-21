


def show_plots():
    import matplotlib.pyplot as plt
    plt.show()


# plot C-P lambda curve
def plotCPowerTsr(sol, show_now = False, add_to_last_plot = False):
    import matplotlib.pyplot as plt
    import numpy as np
    if not add_to_last_plot:
        plt.figure(num="Power-Lambda")
    TSR = []
    CP = []
    
    for i in sol:
        for j in i:
            TSR.append(j.tip_speed_ratio)
            CP.append(j.power_coefficient)
 
    arrTSR = np.array(TSR)
    arrCP = np.array(CP)
    p = arrTSR.argsort()
    
    TSR_sorted = arrTSR[p]
    CP_sorted = arrCP[p]
    
    plt.plot(TSR_sorted, CP_sorted, marker = "4", mec = "black")
    plt.hlines([16/27],linestyles="--", xmin=0, xmax=max(TSR)
               ,colors="red", label="Betz")
    rng = max(TSR) - min(TSR)
    plt.xlim([min(arrTSR)-0.1*rng, max(arrTSR)+0.1*rng])
    plt.grid(True)
    plt.xlabel("Tip Speed Ratio")
    plt.ylabel("Power Coefficient")
    plt.title(r"$C_P$ Vs $\lambda$")
    plt.legend()
    
    
    if show_now == True:
        plt.show()
    
    
def plotCTorqueTsr(sol, show_now = False, add_to_last_plot = False):
    import matplotlib.pyplot as plt
    import numpy as np
    if not add_to_last_plot:
        plt.figure(num="Torque-Lambda")
    TSR = []
    CP = []
    
    for i in sol:
        for j in i:
            TSR.append(j.tip_speed_ratio)
            CP.append(j.torque_coefficient)
    
    arrTSR = np.array(TSR)
    arrCP = np.array(CP)
    p = arrTSR.argsort()
    
    TSR_sorted = arrTSR[p]
    CP_sorted = arrCP[p]
    
    plt.plot(TSR_sorted, CP_sorted, marker = "4", mec = "black")
    plt.xlim([min(arrTSR)-0.1, max(arrTSR)+0.1])
    plt.grid()
    plt.xlabel("Tip Speed Ratio")
    plt.ylabel("Toruqe Coefficient")
    plt.title(r"$C_Q$ Vs $\lambda$")
    plt.legend()
    
    if show_now == True:
        plt.show()
    
    
    
# plot C-P lambda curve
def plotCThrustTsr(sol, show_now = False, add_to_last_plot = False):
    import matplotlib.pyplot as plt
    import numpy as np
    if not add_to_last_plot:
        plt.figure(num="Thrust-Lambda")
    TSR = []
    CT = []
    
    for i in sol:
        for j in i:
            TSR.append(j.tip_speed_ratio)
            CT.append(j.thrust_coefficient)
    
    arrTSR = np.array(TSR)
    arrCT = np.array(CT)
    p = arrTSR.argsort()
    
    TSR_sorted = arrTSR[p]
    CT_sorted = arrCT[p]
    
    plt.plot(TSR_sorted, CT_sorted, marker = "4", mec = "black")
    plt.xlim([min(arrTSR)-0.1, max(arrTSR)+0.1])
    plt.grid(True)
    plt.hlines([8/9],linestyles="--", xmin=0, xmax=max(TSR)
               ,colors="red", label="Betz")
    plt.xlabel("Tip Speed Ratio")
    plt.ylabel("Thrust Coefficient")
    plt.title(r"$C_T$ Vs $\lambda$")
    plt.legend()
    
    if show_now == True:
        plt.show()
    
    
    
def plotAllCoeffsTsr(sol, show_now = False):    
    import matplotlib.pyplot as plt
    import numpy as np
    fig, (ax1, ax2, ax3) = plt.subplot(1,3, num="All Coefficients-Lambda")
    
    # FINISH THIS!
    TSR = []
    CT = []
    CP = []
    CQ = []
    
    for i in sol:
        for j in i:
            TSR.append(j.tip_speed_ratio)
            CT.append(j.thrust_coefficient)
            CQ.append(j.torque_coefficient)
            CP.append(j.power_coefficients)
    
    arrTSR = np.array(TSR)
    arrCT = np.array(CT)
    arrCQ = np.array(CQ)
    arrCP = np.array(CP)
    
    p = arrTSR.argsort()
    
    TSR_sorted = arrTSR[p]
    CT_sorted = arrCT[p]
    CQ_sorted = arrCQ[p]
    CP_sorted = arrCP[p]
    
    plt.plot(TSR_sorted, CT_sorted, marker = "4", mec = "black")
    plt.xlim([min(arrTSR)-0.1, max(arrTSR)+0.1])
    plt.grid()
    plt.hlines([8/9],linestyles="--", xmin=0, xmax=max(TSR)
               ,colors="red", label="Betz")
    plt.xlabel("Tip Speed Ratio")
    plt.ylabel("Thrust Coefficient")
    plt.title("C_T Vs TSR")
    plt.legend()
    plt.grid(True)
    
    if show_now == True:
        plt.show()
    
        

    
def plotBladeLoads(result, show_now = False):  
    import matplotlib.pyplot as plt
    import numpy as np  
    fig, axs = plt.subplots(ncols=1, nrows=3)
    
    radial_positions = result.sectional_positions
    chord_lengths = result.sectional_chord_lengths
    
    thrust = result.sectional_thrust
    torque = result.sectional_torque
    
    total_thrust = result.thrust
    total_torque = result.torque
    total_power = result.power
    
    axs[0].plot(radial_positions, chord_lengths,marker="+")
    axs[0].fill_between(radial_positions,chord_lengths, hatch="x", alpha=0.3)
    axs[0].grid()
    axs[0].set_ylim([-0.5, max(chord_lengths)+1])
    axs[0].set_title("Blade Chord")
    
    axs[1].plot(radial_positions, thrust, marker="+")
    axs[1].grid()
    axs[1].set_title(f"Thrust [Net {total_thrust:.5}N]")

    
    axs[2].plot(radial_positions, torque, marker="+")
    axs[2].grid()
    axs[2].set_title(f"Torque Net [{total_torque:.5}Nm]")
    axs[2].set_xlabel(r"radius $m$")
    
    fig.suptitle(rf"$U_\infty$:{result.wind_speed:<5.5}ms⁻¹     $\Omega$:{result.rot_speed:<5.5}rads⁻¹")
    
    
    plt.subplots_adjust(hspace=0.5)

    if show_now:
        plt.show()
    
    
    
    
    
    
    
    
   