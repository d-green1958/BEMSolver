import source
import simulation
import gui
import case
import plotter


def bempy_path():
    from os import path
    return path.abspath(__path__[0])

def reference_case_path():
    from os import path, listdir, walk, getcwd
    folder_name = "reference_cases"
    folder_path = path.abspath(path.join(__path__[0], folder_name))
    return folder_path