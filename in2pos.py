from os import path
from kit.vasp import POSCAR
from kit.interface import args

if __name__ == "__main__":
    while(1):
        kind = input("Software, QE or Conquest (q/c): ")
        if kind.lower() == "q" or kind.lower() == "c":
            kind = kind.lower()
            break
        else:
            print("Warning: Input 'q' or 'c'.")
    if kind == 'q':
        from kit.software import QE_Single_Point
        arg = args({"output": "POSCAR"}, "QE.in")
        qe = QE_Single_Point(arg)
        poscar = POSCAR(arg)
        qe.read_all()
        poscar.bridge = qe
        poscar.unit_conversion = qe.unit_conversion
        poscar.write_all()
    elif kind == 'c':
        from kit.software import Conquest_Single_Point
        arg = args({"output": "POSCAR", "elementlist": "Conquest_input"}, "ionpos.dat")
        if not path.isfile(arg.args.elementlist):
            if path.isfile(path.join(path.dirname(arg.args.input), "Conquest_input")):
                arg.args.elementlist = path.join(path.dirname(arg.args.input), "Conquest_input")
            elif path.isfile(path.join(arg.args.input, "Conquest_input")):
                arg.args.elementlist = path.join(arg.args.input, "Conquest_input")
        arg.args.elementlist = arg.input_file_check(arg.args.elementlist, "Conquest_input", mandatory=True, isdir=False)
        conquest = Conquest_Single_Point(arg)
        poscar = POSCAR(arg)
        conquest.read_all()
        poscar.bridge = conquest
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
    print("Done!")
