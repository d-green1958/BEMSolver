import source
import simulation
import gui
import case
import plotter
import exceptions
import reference_cases


def bempy_path():
    from os import path
    return path.abspath(__path__[0])