from os import getcwd
from kit.vasp import XDATCAR
from kit.interface import args, single, wrap, smooth
from kit.etc import write_csv

if __name__ == "__main__":
    arg = args({"output": None, "atomfile2": [None, str, "provide the second atom information. default is None"]}, "XDATCAR")
    xdatcar = XDATCAR(arg)
    
    atom_list = single(xdatcar, 1, name="Atom 1")
    xdatcar.args.atomfile = arg.args.atomfile2
    atom_list.put(single(xdatcar, 1, name="Atom 2"))
    
    xdatcar.args = arg.args
    
    xdatcar.read_all(atoms=atom_list)
    
    dist = wrap(xdatcar, type="d")
    
    line, factor = smooth(dist)
    atom_list = atom_list.get(True)
    cen_atom, mea_atom = atom_list[0], atom_list[1]
    filename = (f"{cen_atom} {mea_atom}" if arg.args.output is None else arg.args.output)
    filename = arg.same_name(getcwd(), f"{filename}.csv", factor)
    write_csv(filename, line, tuple(range(xdatcar.final_step)), f"{cen_atom}_{mea_atom}")
    print("Done!")
