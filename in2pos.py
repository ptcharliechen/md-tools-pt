from kit.args import Args
arg = Args({"input": "", "output": "POSCAR", "mode": ["None", str, "the input file mode. c: CONQUEST, g: Gromacs, l: LAMMPS, q: Quantum ESPRESSO"]}, steps=False, atoms=False)
from os import path, getcwd, listdir
from kit.vasp import POSCAR

if __name__ == "__main__":
    arg.same_name(getcwd(), "POSCAR")
    if arg.args.input is None:
        arg.args.input = ""
    if arg.args.mode.lower() not in ["c", "g", "l", "q"] or (not path.isfile(arg.args.input) and not path.isdir(arg.args.input)):
        for f in listdir(getcwd()):
            if "ionpos.dat" in f:
                arg.args.input = f
                if arg.args.elements is None and path.isfile("Conquest_input"):
                    arg.args.elements = "Conquest_input"
                break
            elif "data.lmp" in f:
                arg.args.input = f
                if arg.args.elements is None and path.isfile("in.lmp"):
                    arg.args.elements = "in.lmp"
                break
            elif f == "QE.in" or "gro" in f:
                arg.args.input = f
                break
        else:
            arg.args.input = arg.input_file_check(getcwd(), arg.args.input)
    elif path.isfile(arg.args.input):
        arg.args.input = arg.args.input
    elif path.isdir(arg.args.input):
        for f in listdir(arg.args.input):
            if f == "ionpos.dat":
                arg.args.input = path.join(arg.args.input, f)
                if arg.args.elements is None and path.isfile(arg.args.input, "Conquest_input"):
                    arg.args.elements = path.join(arg.args.input, "Conquest_input")
                break
            elif f == "QE.in":
                arg.args.input = path.join(arg.args.input, f)
                break
            elif f == "data.lmp":
                arg.args.input = path.join(arg.args.input, f)
                if arg.args.elements is None and path.isfile("in.lmp"):
                    arg.args.elements = "in.lmp"
                break
            elif "gro" in f:
                arg.args.input = path.join(arg.args.input, f)
                break
        else:
            arg.args.input = arg.input_file_check(arg.args.input, "")
    else:
        arg.args.input = arg.input_file_check(getcwd(), arg.args.input)
    
    while(1):
        if "QE.in" in arg.args.input:
            arg.args.mode = 'q'
            break
        elif "ionpos.dat" in arg.args.input:
            arg.args.mode = 'c'
            break
        elif "gro" in arg.args.input:
            arg.args.mode = 'g'
            break
        elif "lmp" in arg.args.input:
            arg.args.mode = 'l'
            break
        else:
            kind = input("Iutput format. CONQUEST (c), Gromacs (g), LAMMPS (l), or QE (q): ")
            if kind.lower() == 'q':
                arg.input_type = "QE.in"
            elif kind.lower() == 'c':
                arg.input_type = "ionpos.dat"
            elif kind.lower() == 'g':
                arg.input_type = "gro"
            elif kind.lower() == 'l':
                arg.input_type = "data.lmp"
            else:
                print("Warning: Input 'c', 'g', 'l', or 'q'.")
                continue
            arg.args.mode = kind.lower()
            arg.input_check()
            break
    
    if arg.args.mode == 'q':
        from kit.software import QE_Single_Point
        arg.args.input = arg.input_file_check(arg.args.input, "QE.in")
        qe = QE_Single_Point(arg)
        poscar = POSCAR(arg)
        qe.read_all()
        poscar.bridge = qe
        poscar.unit_conversion = qe.unit_conversion
        poscar.write_all()
    elif arg.args.mode == 'c':
        from kit.software import Conquest_Single_Point
        arg.args.input = arg.input_file_check(arg.args.input, "ionpos.dat")
        if arg.args.elements is None:
            arg.args.elements = arg.input_file_check(getcwd(), "Conquest_input")
        conquest = Conquest_Single_Point(arg)
        conquest.read_all()
        poscar = POSCAR(conquest)
        while(1):
            unit = input("Length unit of source file, angstrom or bohr (a/b) [a]: ")
            if unit == "":
                unit = 'a'
            if unit.lower() == 'b':
                poscar.unit_conversion = 0.529177249
                break
            elif unit.lower() == 'a':
                break
            else:
                print("Warning: Input 'a' or 'b'.")
        poscar.write_all()
    elif arg.args.mode == 'l':
        from kit.software import LAMMPS_Single_Point
        arg.args.input = arg.input_file_check(arg.args.input, "data.lmp")
        if arg.args.elements is None:
            arg.args.elements = arg.input_file_check(getcwd(), "in.lmp")
        lammps = LAMMPS_Single_Point(arg)
        lammps.read_all()
        poscar = POSCAR(lammps)
        poscar.write_all()
    elif arg.args.mode == 'g':
        from kit.software import Gromacs_Single_Point
        arg.args.input = arg.input_file_check(arg.args.input, "gro")
        if arg.args.elements is None:
            arg.args.elements = arg.input_file_check(getcwd(), "topol.top")
        gromacs = Gromacs_Single_Point(arg)
        gromacs.read_all()
        poscar = POSCAR(gromacs)
        poscar.unit_conversion = 10
        poscar.write_all()
    print("Done!")
