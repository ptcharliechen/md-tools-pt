from os import path, mkdir, getcwd
from shutil import copy
from copy import deepcopy
from numpy import array
from kit.vasp import XDATCAR, POSCAR
from kit.interface import step, args

if __name__ == "__main__":
    arg = args({"output": "POSCARs", "incar": [None, str, "INCAR path"], "potcar": [None, str, "POTCAR path"], "kpoints": [None, str, "KPOINTS path"]}, "XDATCAR")
    
    xdatcar = XDATCAR(arg)
    steps = step(xdatcar)
    xdatcar.cutoff_step = max(steps.get())
    xdatcar.read_fast(steps=steps)
    
    if arg.args.incar is None or arg.args.kpoints is None or arg.args.potcar is None:
        while(1):
            yn = input("Would you like to import INCAR, KPOINTS, or POTCAR file (y/n) [n]: ")
            if yn == "":
                yn = 'n'
            if yn.lower() == 'y':
                for exceptPOSCARfile in ["incar", "kpoints", "potcar"]:
                    setattr(arg.args, exceptPOSCARfile, arg.input_file_check(getattr(arg.args, exceptPOSCARfile), exceptPOSCARfile.upper(), mandatory=False, isdir=False))
                break
            elif yn.lower() == 'n':
                for exceptPOSCARfile in ["incar", "kpoints", "potcar"]:
                    if getattr(arg.args, exceptPOSCARfile) is None:
                        print("Warning: Path of {} does not exist.".format(exceptPOSCARfile.upper()))
                break
            else:
                print("Warning: Input 'y' or 'n'.")
    
    mkdir(path.join(getcwd(), arg.args.output))
    for idx, step in enumerate(steps.get(True)):
        mkdir(path.join(getcwd(), arg.args.output, str(step)))
        poscar = POSCAR(arg)
        poscar_arg = deepcopy(arg.args)
        poscar_arg.output = path.join(getcwd(), arg.args.output, str(step), "POSCAR")
        poscar.args, poscar.title, poscar.elements, poscar.atoms = poscar_arg, "POSCAR_{}".format(idx), xdatcar.elements, xdatcar.atoms
        poscar.lattice = (xdatcar.lattice[step-arg.arg.step] if array(xdatcar.lattice).ndim == 3 else xdatcar.lattice)
        poscar.fast_position = xdatcar.fast_position[idx]
        poscar.write_fast()
        for exceptPOSCARfile in ["incar", "kpoints", "potcar"]:
            if getattr(arg.args, exceptPOSCARfile) is not None:
                copy(getattr(arg.args, exceptPOSCARfile), path.join(getcwd(), arg.args.output, str(step)))
    print("Done!")
