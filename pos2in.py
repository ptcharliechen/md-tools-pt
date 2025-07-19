from kit.vasp import POSCAR
from kit.interface import args

if __name__ == "__main__":
    while(1):
        kind = input("Output software, QE or Conquest (q/c): ")
        if kind.lower() == "q" or kind.lower() == "c":
            kind = kind.lower()
            break
        else:
            print("Warning: Input 'q' or 'c'.")
    if kind == 'q':
        from kit.software import QE_Single_Point
        arg = args({"output": "QE.in"}, "POSCAR")
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
    elif kind == 'c':
        from kit.software import Conquest_Single_Point
        arg = args({"output": "ionpos.dat"}, "POSCAR")
        poscar = POSCAR(arg)
        conquest = Conquest_Single_Point(arg)
        poscar.read_all()
        conquest.bridge = poscar
        conquest.write_all()
    print("Done!")
