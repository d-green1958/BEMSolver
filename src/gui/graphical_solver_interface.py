# from ..BEMsolver import *
import tkinter as tk
from tkinter import messagebox
from tkinter import ttk
import customtkinter as ctk

def create_tab(notebook, text):
    frame = ttk.Frame(notebook)
    notebook.add(frame, text=text)
    return frame

ctk.set_appearance_mode("System")  
ctk.set_default_color_theme("green") 

class GraphicalSolver(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("BEM Solver")    
        
        
        self.options_frame = ctk.CTkFrame(self, width=400, corner_radius=30)
        self.options_frame.grid(row=0, column=1, sticky="nsew", padx=20, pady = 20)
        self.options_label = ctk.CTkLabel(self.options_frame, text="Options", font=ctk.CTkFont(size=20, weight="bold"))
        self.options_label.grid(row=0, column=0, padx=20, pady = 20)
        
        self.tip_loss_model_box = ctk.CTkComboBox(self.options_frame, values=["none", "Prandtl approx", "Prandtl", "Xu & Sankar"])
        self.tip_loss_model_box.grid(row=1, column=1, padx=20, pady = 20)
        self.tip_loss_model_label = ctk.CTkLabel(self.options_frame, text="Tip loss model")
        self.tip_loss_model_label.grid(row=1, column=0, pady=10, padx=10)
        
        self.hub_loss_model_box = ctk.CTkComboBox(self.options_frame, values=["none", "Prandtl approx", "Prandtl"])
        self.hub_loss_model_label = ctk.CTkLabel(self.options_frame, text="Hub loss model")
        self.hub_loss_model_label.grid(row=2, column=0, pady=10, padx=10)
        self.hub_loss_model_box.grid(row=2, column=1, pady=10, padx=10)     
        
        self.results_frame = ctk.CTkFrame(self, width=400, corner_radius=30)
        self.results_frame.grid(row=1, column=0, columnspan=2, sticky="nsew", padx=20, pady = 20)
        self.results_label = ctk.CTkLabel(self.results_frame, text="Results", font=ctk.CTkFont(size=20, weight="bold"))
        self.results_label.grid(row=0, column=0, padx=20, pady = 20)
        
        

        self.tabview = ctk.CTkTabview(self, width=400)
        self.tabview.grid(row=0, column=0, padx=10, pady = 20, sticky="nsew")
           
        self.tabview.add("Single")
        self.tabview.add("Parametric")
        
        # SINGLE RUN TAB
        self.wind_speed_label_single = ctk.CTkLabel(self.tabview.tab("Single"), text="Wind Speed [ms^-1]")
        self.wind_speed_label_single.grid(row=0, column=0, pady = 20, padx = 20)
        self.wind_speed_single = ctk.CTkEntry(self.tabview.tab("Single"),)
        self.wind_speed_single.grid(row = 0, column = 1, pady = 20, padx = 20)
        
        self.rot_speed_label_single = ctk.CTkLabel(self.tabview.tab("Single"), text="Rotational Speed [rads^-1]")
        self.rot_speed_label_single.grid(row=1, column=0, pady = 20, padx = 20)
        self.rot_speed_single = ctk.CTkEntry(self.tabview.tab("Single"),)
        self.rot_speed_single.grid(row = 1, column = 1, pady = 20, padx = 20)
    
        self.run_button_single = ctk.CTkButton(self.tabview.tab("Single"), text="Run", command=self.run_single)       
        self.run_button_single.grid(row=5, column=2)
        # self.run_button_single.pack()
        
        # PARAMETRIC RUN TAB
        
    def run_single(self):
        wind_speed = self.wind_speed_single.get()
        rot_speed = self.rot_speed_single.get()
        if not wind_speed.isdigit():
            print("Err: invalid wind speed")
            messagebox.showerror("Invalid Entry")
            raise TypeError
        else:
            wind_speed = float(wind_speed)
            
        if not rot_speed.isdigit():
            print("Err: invalid rotational speed")
            messagebox.showerror("Invalid Entry")
            raise TypeError
        else:
            rot_speed = float(rot_speed)
        
        if (wind_speed < 0) or (rot_speed < 0):
            messagebox.showerror("Invalid Entry")
            raise ValueError
        
        
        return
        
        
        
        
        


if __name__ == "__main__":
    app = GraphicalSolver()
    app.mainloop()

