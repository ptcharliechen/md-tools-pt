from os import path, getcwd, listdir
from numpy import sum, zeros, linspace
from scipy.interpolate import interp1d
from kit.fundamental import Step
from kit.vasp import ACF, ACFs, Charge_Edit, Args
from kit.interface import lines, smooth
from kit.etc import write_csv

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
        modeFlag = False
        if not modeFlag:
            mode = input("Single point or trajectory mode (s/t): ")
        args.args.input = path.abspath(args.args.input)
        if mode.lower() == 's':
            modeFlag = True
            args.input_type = "ACF.dat"
            args.input_check()
            if path.isfile(args.args.input):
                return args
            elif path.isdir(args.args.input):
                return args
        elif mode.lower() == 't':
            modeFlag = True
            flag = False
            args.args.input = input("Input directory path with ACF.dat files: ")
            for filename in listdir(args.args.input):
                if not filename.isdigit() and path.isfile(path.join(args.args.input, filename, "ACF.dat")):
                    flag = True
                elif not filename.isdigit() and not path.isfile(path.join(args.args.input, filename, "ACF.dat")):
                    input(f"Warning: {path.join(args.args.input, filename)} does not have any ACF.dat file. Press Enter button to continue after checking.")
            if flag:
                return args
        else:
            print("Warning: Input error.")

if __name__ == "__main__":
    arg = Args({"output": "Selected_Charge_Diff", "molecules": [None, str, "specify a file with molecule information"],
                "poscar": [None, str, "POSCAR path to get element information"], "potcar": [None, str, "POTCAR path to get charge information"]})
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
    
    if arg.args.molecules is None:
        print("\nWarning: You do not specify a file with molecule information.")
    elif not path.isfile(arg.args.molecules):
        print(f"\nWarning: {arg.args.molecules} does not exist.")
    if arg.args.molecules is None or not path.isfile(arg.args.molecules):
        print("It is optional. You can choose whether to provide the molecule information or not.")
        while(1):
            tmp = input("\nDo you want to provide the molecule information? (y/n) [n]: ")
            if tmp == "":
                tmp = "n"
            if tmp.lower() == 'y':
                charge_edit = Charge_Edit(arg, steps=Step_Charge(arg))
                charge_edit.read_molecules()
                break
            elif tmp.lower() == 'n':
                charge_edit = Charge_Edit(arg, steps=Step_Charge(arg))
                break
            else:
                print("Warning: Input error.")
    else:
        charge_edit = Charge_Edit(arg, steps=Step_Charge(arg))
        
    while(1):
        try:
            charge_edit.read_elements()
        except:
            arg.args.poscar = input("Input POSCAR file: ")
            charge_edit.args = arg.args
        else:
            break
        charge_edit.read_elements()
    while(1):
        try:
            charge_edit.read_ref()
        except:
            arg.args.potcar = input("Input POTCAR file: ")
            charge_edit.args = arg.args
        else:
            break
        charge_edit.read_elements()
    
    if path.isfile(arg.args.input):
        atom_lists = lines(charge_edit)
        data = []
        
        for atom_list in atom_lists:
            tmp = []
            charge = ACF(arg)
            charge.read_atom(path.join(path.dirname(arg.args.input), "ACF.dat"), atoms=atom_list)
            tmp.append(f"{sum(charge.charge):.4f}")
            tmp.append(charge.charge)
            tmp.append(charge.atoms.get())
            data.append(tmp)
        write_csv(arg.args.output, data, list(range(len(data))), ["charge sum", "charge", "atoms"])
    elif path.isdir(arg.args.input):
        if not path.isfile(path.join(arg.args.input, "Charge.csv")) or not path.isfile(path.join(arg.args.input, "Charge_Diff.csv")):
            charge_edit.build()
            charge_edit.summary("Charge", "Charge_Diff")
        
        atom_lists = lines(charge_edit)
        charges = []
        for atom_list in atom_lists:
            charge = ACFs(charge_edit, steps=Step_Charge(charge_edit))
            charge.read_atom(atoms=atom_list, diff_file="Charge_Diff")
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

