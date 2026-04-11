from kit.args import Args
arg = Args({"output": "Selected_Charge_Diff", "molecules": ["molecules.csv", str, "specify a file with molecule information"], \
            "poscar": [None, str, "POSCAR path to get element information"], "potcar": [None, str, "POTCAR path to get electron information"]})
from os import path, getcwd, listdir
from numpy import zeros, linspace
from scipy.interpolate import interp1d
from kit.fundamental import Step
from kit.vasp import ACF, ACFs, Charge_Edit
from kit.interface import lines, smooth

class Step_Charge(Step):
    def __init__(self, input_obj=None, put_flag=False):
        Step.__init__(self, input_obj, put_flag)
    @property
    def args(self):
        return self._args
    @args.setter
    def args(self, args):
        self._args = args
        self._start = self._args.step

def choose_mode(args):
    while(1):
        mode_flag = False
        mode = input("Single point or trajectory mode (s/t): ")
        args.args.input = path.abspath(args.args.input) if args.args.input is not None else getcwd()
        if mode.lower() == 's':
            mode_flag = True
            args.input_type = "ACF.dat"
            args.input_check()
            if path.isfile(args.args.input):
                return args
            elif path.isdir(args.args.input):
                return args
        elif mode.lower() == 't':
            mode_flag = True
            flag = False
            dirpath = input("Input directory path with ACF.dat files: ")
            if not path.isdir(dirpath):
                print(f"Warning: {dirpath} is not a directory.")
                continue
            else:
                args.args.input = dirpath
            for filename in listdir(args.args.input):
                if path.isfile(path.join(args.args.input, filename)):
                    continue
                elif filename.isdigit():
                    if path.isfile(path.join(args.args.input, filename, "ACF.dat")):
                        flag = True
                    else:
                        input(f"Warning: {path.join(args.args.input, filename)} does not have any ACF.dat file. Press Enter button to continue after checking.")
            if flag:
                return args
        else:
            print("Warning: Input error.")

if __name__ == "__main__":
    arg.same_name(getcwd(), arg.args.output)

    if arg.args.input is None:
        arg = choose_mode(arg)
    elif path.isfile(arg.args.input) or path.isfile(path.join(arg.args.input, "ACF.dat")):
        arg.args.input = arg.input_file_check(arg.args.input, "ACF.dat")
    else:
        flag = False
        if path.isdir(arg.args.input):
            for filename in listdir(arg.args.input):
                if filename.isdigit() and path.isfile(path.join(arg.args.input, filename, "ACF.dat")):
                    flag = True
                elif filename.isdigit() and not path.isfile(path.join(arg.args.input, filename, "ACF.dat")):
                    input(f"Warning: {path.join(arg.args.input, filename)} does not have any ACF.dat file. Press Enter button to continue after checking.")
        if not flag:
            arg = choose_mode(arg)
    
    charge_edit = Charge_Edit(arg, steps=Step_Charge(arg))
    if arg.args.molecules is not None:
        charge_edit.read_molecules()
    
    while(1):
        try:
            charge_edit.read_elements()
        except:
            filepath = input("Input POSCAR file: ")
            if path.isfile(path.join(filepath, "POSCAR")):
                filepath = path.join(filepath, "POSCAR")
            if not path.isfile(filepath):
                print(f"Warning: {filepath} does not exist.")
            else:
                arg.args.poscar = filepath
                charge_edit.args = arg.args
                charge_edit.read_elements()
                break
        else:
            break
    
    while(1):
        try:
            charge_edit.read_ref()
        except:
            filepath = input("Input POTCAR file: ")
            if path.isfile(path.join(filepath, "POTCAR")):
                filepath = path.join(filepath, "POTCAR")
            if not path.isfile(filepath):
                print(f"Warning: {filepath} does not exist.")
            else:
                arg.args.potcar = filepath
                charge_edit.args = arg.args
                charge_edit.read_ref()
                break
        else:
            break
    
    if path.isfile(arg.args.input):
        atom_lists = lines(charge_edit)
        data = []
        for atom_list in atom_lists:
            tmp = []
            charge = ACF(charge_edit)
            charge.read_all(atoms=atom_list)
            charge.write_all()
    elif path.isdir(arg.args.input):
        if not path.isfile(path.join(arg.args.input, "Charge.csv")) or not path.isfile(path.join(arg.args.input, "Charge_Diff.csv")):
            charge_edit.build()
            charge_edit.summary("Charge", "Charge_Diff")
        
        atom_lists = lines(charge_edit)
        charges = []
        for atom_list in atom_lists:
            charge = ACFs(charge_edit, steps=Step_Charge(charge_edit))
            charge.read_all(atoms=atom_list, diff_file="Charge_Diff")
            charge.write_all()
            charges.append(charge)
        
        steps = charge.steps.get()
        lines = zeros((max(steps), len(charges)))
        for idx, charge in enumerate(charges):
            fun = interp1d(steps, charge.diff.sum(axis=1), kind="cubic")
            lines[:, idx] = fun(linspace(steps[0], steps[-1], lines.shape[0]))
        
        smoothed_lines, smoothed_factor = smooth(lines)
        if not smoothed_factor:
            filename = f"{arg.args.output}_Continuity.csv"
        else:
            filename = f"{arg.args.output}_Alpha={smoothed_factor}.csv"
        arg.same_name(getcwd(), filename)
        with open(filename, "w") as write_file:
            for step in range(arg.args.step, max(steps)+arg.args.step):
                write_file.write("{}".format(step))
                for _ in smoothed_lines[step-arg.args.step]:
                    write_file.write(",{:.4f}".format(_))
                write_file.write('\n')
    print("Done!")
