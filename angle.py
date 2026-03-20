from kit.args import args
arg = args({"output": "", "radian": [False, "use radian as output unit"]}, "XDATCAR", elements=False)
from os import getcwd
from kit.vasp import XDATCAR
from kit.interface import single, wrap, smooth
from kit.etc import write_csv

if __name__ == "__main__":
    xdatcar = XDATCAR(arg)
    
    atoms = single(xdatcar, 1, line=[1], name="Center atom")
    atoms.put(single(xdatcar, 2, line=[2], show_info_flag=False, name="Side atoms"))
    
    xdatcar.read_all(atoms=atoms)
    angle = wrap(xdatcar, type="a")
    
    line, factor = smooth(angle)
    
    atoms = atoms.get(True)
    cen_atom, side_atom_1, side_atom_2 = atoms[0], atoms[1], atoms[2]
    
    filename = (f"{cen_atom}_{side_atom_1}_{side_atom_2}" if arg.args.output == "" else arg.args.output)
    filename = arg.same_name(getcwd(), filename+".csv", factor)
    write_csv(filename, line, tuple(range(xdatcar.final_step)), f"{cen_atom}_{side_atom_1}_{side_atom_2}")
    print("Done!")
