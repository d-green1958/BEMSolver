import os

def make_case():
    cDir = os.getcwd()
    
    list = os.listdir()

    if "turbineProps" not in list:
        print(">> building turbineProps dir")
        os.mkdir("turbineProps")
        
        print(">> building aerofoilData dir")
        os.mkdir("turbineProps/aerofoilData")
        
        
    if "simProps" not in list:
        print(">> building simProps dir")
        os.mkdir("simProps")
    
    if "results" not in list:
        print(">> building results dir")
        os.mkdir("results")
        
    ## FINISH THIS!!!
        
