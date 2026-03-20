from os import path
from collections.abc import Iterable
from argparse import ArgumentParser, Namespace
from copy import deepcopy

class Args:
    def __init__(self, default=None, arg_func=None, steps=True, atoms=True, elements=True):
        args = ArgumentParser()
        if arg_func is None:
            arg_func = self.func
        
        args = arg_func(default, args, steps, atoms, elements)
        
        self.__args = args.parse_args()
    def func(self, default, args, steps, atoms, elements):
        args.add_argument("input", type=str, nargs='?', help="input file or directory path")
        if steps:
            args.add_argument("-step", dest="step", type=int, default=1, help="initial step number of the trajectory file. default is 1.")
            args.add_argument("-stepfile", dest="stepfile", type=str, default=None, help="specify a file with the step information. default is None.")
        if atoms:
            args.add_argument("-atom", dest="atom", type=int, default=0, help="initial atom number. default is 0.")
            args.add_argument("-atomfile", dest="atomfile", type=str, default=None, help="specify a file with the atom information. default is None.")
        if elements:
            args.add_argument("-elements", dest="elements", type=str, default=None, help="specify a file with the element information. default is None.")
            args.add_argument("-elementfile", dest="elementfile", type=str, default=None, help="specify a file with self-defined elements name. default is None.")
        if isinstance(default, dict):
            for key in default.keys():
                if key == "output" and not isinstance(default[key], list):
                    args.add_argument("-output", dest="output", type=str, default=default[key], help=f"output file or directory name. default is {default[key]}.")
                elif isinstance(default[key], Iterable) and len(default[key]) == 3:
                    args.add_argument(f"-{key}", dest=key, type=default[key][1], default=default[key][0], help=f"{default[key][2]}. default is {default[key][0]}.")
                elif isinstance(default[key], Iterable) and len(default[key]) == 2:
                    args.add_argument(f"-{key}", dest=key, type=default[key][0].__class__, default=default[key][0], help=f"{default[key][1]}. default is {default[key][0]}.")
                elif isinstance(default[key], Iterable) and isinstance(default[key], str):
                    args.add_argument(f"-{key}", dest=key, type=default[key].__class__, default=default[key], help=f"{default[key]}")
                elif isinstance(default[key], Iterable):
                    args.add_argument(f"-{key}", dest=key, type=default[key][0].__class__, default=default[key][0], help=f"{default[key][0]}")
                else:
                    args.add_argument(f"-{key}", dest=key, type=default[key].__class__, default=default[key], help=f"{key}")
        elif default is not None:
            raise ValueError("default should be a dictionary.")
        return args
    @property
    def args(self):
        return self.__args
    @property
    def input_type(self):
        return self.__input_type
    @input_type.setter
    def input_type(self, input_type):
        self.__input_type = input_type
    @staticmethod
    def same_name(present_path, filename, factor=0, isdir=False):
        while(1):
            if "_Alpha=" in filename:
                repeat = True
            else:
                repeat = False
            if filename is None:
                break
            if factor and not repeat:
                tmp = path.splitext(filename)
                filename = tmp[0] + f"_Alpha={factor}" + tmp[1]
            file = path.join(present_path, filename)
            if not isdir and path.isfile(file):
                action = input(f"Warning: '{filename}' file exists in '{present_path}'. Delete, move, or rename (d/m/n): ")
                repeat = True
            elif not isdir and path.isfile(f"{file}.csv"):
                action = input(f"Warning: '{filename}.csv' file exists in '{present_path}'. Delete, move, or rename (d/m/n): ")
                filename += ".csv"
                repeat = True
            elif not isdir and path.isfile(f"{file}.dat"):
                action = input(f"Warning: '{filename}.dat' file exists in '{present_path}'. Delete, move, or rename (d/m/n): ")
                filename += ".dat"
                repeat = True
            elif isdir and path.isdir(file):
                action = input(f"Warning: '{filename}' directory exists in '{present_path}'. Delete, move, or rename (d/m/n): ")
                repeat = True
            else:
                return filename
            Args.action(present_path, filename, action)
    def input_check(self):
        self.__args.input = self.input_file_check(self.__args.input, self.__input_type)
    def action(present_path, filename, action):
        file_path = path.join(present_path, filename)
        if action.lower() == 'd' and path.isdir(file_path):
            from shutil import rmtree
            rmtree(file_path)
        elif action.lower() == 'd' and path.isfile(file_path):
            from os import remove
            remove(file_path)
        elif action.lower() == 'n':
            from os import rename
            newName = input("New name: ")
            rename(file_path, path.join(present_path, newName))
        elif action.lower() == 'm':
            from shutil import move
            new_path = input("New path: ")
            if path.isdir(new_path):
                move(file_path, new_path)
            else:
                print("Warning: Provide a path to existence.")
    @staticmethod
    def input_file_check(input_path, input_type, mandatory=True, isdir=False):
        while(1):
            if input_path is None and not mandatory:
                input_path = input(f"Input {input_type} path (If not needed, input 'no'): ")
            elif input_path is None:
                input_path = input(f"Input {input_type} path: ")
            if not mandatory and input_path.lower() == "no":
                return None
            if not isdir:
                if path.isfile(input_path):
                    return input_path
                elif path.isfile(path.join(input_path, input_type)):
                    return path.join(input_path, input_type)
                elif path.isfile(input_path+input_type):
                    return input_path+input_type
                elif path.isfile(input_path+"."+input_type):
                    return input_path+"."+input_type
                else:
                    print(f"Warning: '{input_type}' does not exist in {path.abspath(input_path)}.\n")
                    input_path = None
            else:
                if path.isdir(input_path):
                    return input_path
                elif path.isdir(path.join(input_path, input_type)):
                    return path.join(input_path, input_type)
                else:
                    print(f"Warning: {input_type} does not exist in {path.abspath(input_path)}.")
                    input_path = None
    @staticmethod
    def arg_check(args):
        if isinstance(args, dict):
            args = Args(args)
        if isinstance(args, Args):
            return deepcopy(args.args)
        elif isinstance(args, Namespace):
            return deepcopy(args)
        else:
            return None

def args(parameters: dict, input_type: str, molecules=False, steps=True, atoms=True, elements=True, input_isdir=False, output_isdir=False) -> Args:
    """
    Check and set the input and output files for preprocessing.

    input check:
        The input file exists or not. If not, reset the path.

    Parameters
    ----------
    parameters (mandatory): dict
        A dictionary of the parameter options allowed by the designer.
    input_type (mandatory): str
        Filename extension or default input file name. Choose one.
        'input_file_check' method helps you judge it.
        Attention, default input file name has higher priority.
        If not, filename extension is applied.
        Default input file name:
            If users just provide a directory, 'input_file_check' method implements it
            with 'input_type'.
        Filename extension:
            If users give a file without filename extension, 'input_file_check' method
            supplements it with 'input_type'.

    Returns
    -------
    args : Args
        A dictionary of preprocessing options.
    """
    from os import listdir, getcwd

    # Initialize Args with parameters
    if molecules:
        parameters["molecules"] = [None, str, "specify a file with molecule information"]
        if "molecules.csv" in listdir():
            parameters["molecules"][0] = "molecules.csv"
    args = Args(parameters, steps=steps, atoms=atoms, elements=elements)
    
    if args.args.input is None:
        for file in listdir():
            if input_type in file:
                args.args.input = file
                break

    # Check the 'input' key in 'parameters' dictionary is complete or not.
    # If not, a valid file path should be reset to 'input' key in 'args' by users.
    args.args.input = args.input_file_check(args.args.input, input_type, isdir=input_isdir)
    
    # Check the 'output' key in 'parameters' dictionary exists any same name file or directory.
    # If so, the old file or directory name should be revised by users.
    if not hasattr(args.args, "samename") or args.args.samename:
        args.same_name(getcwd(), args.args.output, isdir=output_isdir)
    
    return args
