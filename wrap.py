from kit.args import args
arg = args({"output": ""}, "XDATCAR", elements=False)
from os import getcwd
from kit.vasp import XDATCAR
from kit.interface import single, wrap, smooth
from kit.etc import write_csv

if __name__ == "__main__":
    xdatcar = XDATCAR(arg)
    
    atom_list = single(xdatcar, 1, line=[1], name="Central atom")
    atom_list.put(single(xdatcar, 1, line=[2], show_info_flag=False, name="Measured atom"))
    
    xdatcar.read_all(atoms=atom_list)
    
    dist = wrap(xdatcar, type="d")

    line, factor = smooth(dist[0] if dist.ndim == 2 else dist[0][0].T)
    line = line.T if line.ndim == 2 else line
    atom_list = atom_list.get(True)
    cen_atom, mea_atom = atom_list[0], atom_list[1]
    filename = (f"{cen_atom}_{mea_atom}" if arg.args.output == "" else arg.args.output)
    filename = arg.same_name(getcwd(), filename+".csv", factor)
    write_csv(filename, line if line.ndim < 3 else line[0][0], tuple(range(xdatcar.final_step)), f"{cen_atom}_{mea_atom}")
    print("Done!")
