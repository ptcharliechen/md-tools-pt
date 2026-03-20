from kit.args import args
arg = args({"output": "None", "mode": ["None", str, "the input file mode. c: CONQUEST, g: Gromacs, l: LAMMPS, q: Quantum ESPRESSO"], \
            "molecules": [None, str, "specify a file with molecule information. Only available on LAMMPS mode and Gromacs mode"]}, "POSCAR", steps=False, elements=False)
from os import getcwd
from kit.vasp import POSCAR

if __name__ == "__main__":
    if arg.args.mode.lower() not in ["c", "g", "l", "q"]:
        while(1):
            kind = input("Output format, CONQUEST (c), Gromacs (g), LAMMPS (l), or QE (q): ")
            if kind.lower() in ["c", "g", "l", "q"]:
                arg.args.mode = kind.lower()
                break
            else:
                print("Warning: Input 'c', 'g', 'l', or 'q'.")
    
    if arg.args.mode == 'q':
        from kit.software import QE_Single_Point
        if arg.args.output == "None":
            arg.args.output = arg.same_name(getcwd(), "QE.in")
        poscar = POSCAR(arg)
        while(1):
            tmp = input("Unit of QE output file, fractional, angstrom, or bohr (f/a/b) [f]: ")
            if tmp == "":
                tmp = "f"
            if tmp.lower() == "f":
                unit = "crystal"
                break
            elif tmp.lower() == "a":
                unit = "angstrom"
                break
            elif tmp.lower() == "b":
                unit = "bohr"
                poscar.unit_conversion = 1/0.529177249
                break
            else:
                print("Warning: Input 'f', 'a', or 'b'.")
        qe = QE_Single_Point(arg)
        qe.unit_conversion = poscar.unit_conversion
        poscar.read_all()
        qe.bridge = poscar
        qe.write_all(unit=unit)
    elif arg.args.mode == 'c':
        from kit.software import Conquest_Single_Point
        if arg.args.output == "None":
            arg.args.output = arg.same_name(getcwd(), "ionpos.dat")
        poscar = POSCAR(arg)
        poscar.read_all()
        while(1):
            tmp = input("Unit of CONQUEST output file, angstrom, or bohr (a/b) [b]: ")
            if tmp == "":
                tmp = "b"
            if tmp.lower() == "a":
                unit = "angstrom"
                break
            elif tmp.lower() == "b":
                unit = "bohr"
                poscar.unit_conversion = 1/0.529177249
                break
            else:
                print("Warning: Input 'a', or 'b'.")
        conquest = Conquest_Single_Point(poscar)
        conquest.unit_conversion = poscar.unit_conversion
        conquest.write_all()
    elif arg.args.mode == 'l':
        from os import path
        from kit.software import LAMMPS_Single_Point
        if arg.args.output == "None":
            arg.args.output = arg.same_name(getcwd(), "data.lmp")
        if arg.args.molecules is None and path.isfile("molecules.csv"):
            arg.args.molecules = "molecules.csv"
        poscar = POSCAR(arg)
        poscar.read_all()
        lammps = LAMMPS_Single_Point(poscar)
        lammps.write_all()
    elif arg.args.mode == 'g':
        from os import path
        from kit.software import Gromacs_Single_Point
        if arg.args.output == "None":
            arg.args.output = arg.same_name(getcwd(), "output.gro")
        if arg.args.molecules is None and path.isfile("molecules.csv"):
            arg.args.molecules = "molecules.csv"
        poscar = POSCAR(arg)
        poscar.read_all()
        gromacs = Gromacs_Single_Point(poscar)
        gromacs.unit_conversion = 0.1
        gromacs.write_all()
    print("Done!")from kit.args import args
arg = args({"output": "None", "mode": ["None", str, "the input file mode. c: CONQUEST, g: Gromacs, l: LAMMPS, q: Quantum ESPRESSO"], \
            "molecules": [None, str, "specify a file with molecule information. Only available on LAMMPS mode and Gromacs mode"]}, "POSCAR", steps=False, elements=False)
from os import getcwd
from kit.vasp import POSCAR

if __name__ == "__main__":
    if arg.args.mode.lower() not in ["c", "g", "l", "q"]:
        while(1):
            kind = input("Output format, CONQUEST (c), Gromacs (g), LAMMPS (l), or QE (q): ")
            if kind.lower() in ["c", "g", "l", "q"]:
                arg.args.mode = kind.lower()
                break
            else:
                print("Warning: Input 'c', 'g', 'l', or 'q'.")
    
    if arg.args.mode == 'q':
        from kit.software import QE_Single_Point
        if arg.args.output == "None":
            arg.args.output = arg.same_name(getcwd(), "QE.in")
        poscar = POSCAR(arg)
        while(1):
            tmp = input("Unit of QE output file, fractional, angstrom, or bohr (f/a/b) [f]: ")
            if tmp == "":
                tmp = "f"
            if tmp.lower() == "f":
                unit = "crystal"
                break
            elif tmp.lower() == "a":
                unit = "angstrom"
                break
            elif tmp.lower() == "b":
                unit = "bohr"
                poscar.unit_conversion = 1/0.529177249
                break
            else:
                print("Warning: Input 'f', 'a', or 'b'.")
        qe = QE_Single_Point(arg)
        qe.unit_conversion = poscar.unit_conversion
        poscar.read_all()
        qe.bridge = poscar
        qe.write_all(unit=unit)
    elif arg.args.mode == 'c':
        from kit.software import Conquest_Single_Point
        if arg.args.output == "None":
            arg.args.output = arg.same_name(getcwd(), "ionpos.dat")
        poscar = POSCAR(arg)
        poscar.read_all()
        while(1):
            tmp = input("Unit of CONQUEST output file, angstrom, or bohr (a/b) [b]: ")
            if tmp == "":
                tmp = "b"
            if tmp.lower() == "a":
                unit = "angstrom"
                break
            elif tmp.lower() == "b":
                unit = "bohr"
                poscar.unit_conversion = 1/0.529177249
                break
            else:
                print("Warning: Input 'a', or 'b'.")
        conquest = Conquest_Single_Point(poscar)
        conquest.unit_conversion = poscar.unit_conversion
        conquest.write_all()
    elif arg.args.mode == 'l':
        from os import path
        from kit.software import LAMMPS_Single_Point
        if arg.args.output == "None":
            arg.args.output = arg.same_name(getcwd(), "data.lmp")
        if arg.args.molecules is None and path.isfile("molecules.csv"):
            arg.args.molecules = "molecules.csv"
        poscar = POSCAR(arg)
        poscar.read_all()
        lammps = LAMMPS_Single_Point(poscar)
        lammps.write_all()
    elif arg.args.mode == 'g':
        from os import path
        from kit.software import Gromacs_Single_Point
        if arg.args.output == "None":
            arg.args.output = arg.same_name(getcwd(), "output.gro")
        if arg.args.molecules is None and path.isfile("molecules.csv"):
            arg.args.molecules = "molecules.csv"
        poscar = POSCAR(arg)
        poscar.read_all()
        gromacs = Gromacs_Single_Point(poscar)
        gromacs.unit_conversion = 0.1
        gromacs.write_all()
    print("Done!")
