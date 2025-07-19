from os import getcwd
from kit.vasp import XDATCAR
from kit.interface import args, single, wrap, smooth
from kit.etc import write_csv

if __name__ == "__main__":
    arg = args({"output": None, "atomfilecollinear": [None, str, "provide the collinear atoms information. default is None"], \
                "atomfileside": [None, str, "provide the side atoms information. default is None"], "radian": False}, "XDATCAR")
    xdatcar = XDATCAR(arg)
    
    atoms = single(xdatcar, 1, "Center atom")
    xdatcar.args.atomfile = arg.args.atomfilecollinear
    atoms.put(single(xdatcar, 1, "Collinear atom"))
    xdatcar.args.atomfile = arg.args.atomfileside
    atoms.put(single(xdatcar, 2, "Side atoms"))
    xdatcar.args = arg.args
    
    xdatcar.read_all(atoms=atoms)
    dist = wrap(xdatcar, type="dh")
    
    line, factor = smooth(dist)
    atoms = atoms.get(True)
    cen_atom, collin_atom, side_atom_1, side_atom_2 = atoms[0], atoms[1], atoms[2], atoms[3]
    
    filename = ("{}_{}_{}_{}".format(cen_atom, collin_atom, side_atom_1, side_atom_2) if arg.args.output is None else arg.args.output)
    arg.same_name(getcwd(), "{}.csv".format(filename), factor)
    
    write_csv(filename, line, tuple(range(xdatcar.final_step)), "{}_{}_{}_{}".format(cen_atom, collin_atom, side_atom_1, side_atom_2))
    print("Done!")
