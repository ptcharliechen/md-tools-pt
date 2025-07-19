from os import getcwd
from kit.vasp import POSCAR
from kit.software import Gaussian
from kit.interface import args

if __name__ == "__main__":
    arg = args({"output": "POSCAR.gjf"}, "gjf")
    poscar = POSCAR(arg)
    poscar.read_all()
    gaussian = Gaussian(arg)
    gaussian.bridge = poscar
    gaussian.args.output = poscar.title
    arg.same_name(getcwd(), "{}.gjf".format(gaussian.args.output))
    gaussian.cartesian_position = poscar.cartesian_position
    gaussian.write_all()
    print("Done!")
