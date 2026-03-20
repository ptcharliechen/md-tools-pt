from kit.args import args
arg = args({"output": "poscar"}, "gjf", steps=False, atoms=False)
from os import getcwd
from kit.vasp import POSCAR
from kit.software import Gaussian

if __name__ == "__main__":
    poscar = POSCAR(arg)
    poscar.read_all()
    gaussian = Gaussian(poscar)
    gaussian.args.output = poscar.title
    arg.same_name(getcwd(), f"{gaussian.args.output}.gjf")
    gaussian.cartesian_position = poscar.cartesian_position
    gaussian.write_all()
    print("Done!")
