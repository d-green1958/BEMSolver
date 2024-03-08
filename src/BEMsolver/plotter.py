import matplotlib.pyplot as plt
import numpy as np

# plot C-P lambda curve
def plotCPowerTsr(sol, show_now = False):
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
    plt.xlim([min(arrTSR)-0.1, max(arrTSR)+0.1])
    plt.grid()
    plt.xlabel("Tip Speed Ratio")
    plt.ylabel("Power Coefficient")
    plt.title("C_P Vs TSR")
    plt.legend()
    
    if show_now == True:
        plt.show()
    
    
def plotCTorqueTsr(sol, show_now = False):
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
    plt.title("C_Q Vs TSR")
    plt.legend()
    
    if show_now == True:
        plt.show()
    
    
    
# plot C-P lambda curve
def plotCThrustTsr(sol, show_now = False):
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
    plt.grid()
    plt.hlines([8/9],linestyles="--", xmin=0, xmax=max(TSR)
               ,colors="red", label="Betz")
    plt.xlabel("Tip Speed Ratio")
    plt.ylabel("Thrust Coefficient")
    plt.title("C_T Vs TSR")
    plt.legend()
    
    if show_now == True:
        plt.show()
    
    
    
    
    
   