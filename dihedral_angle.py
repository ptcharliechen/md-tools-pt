from kit.args import args
arg = args({"output": "", "radian": [False, "use radian as output unit"]}, "XDATCAR", elements=False)
from os import getcwd
from kit.vasp import XDATCAR
from kit.interface import single, wrap, smooth
from kit.etc import write_csv

if __name__ == "__main__":
    xdatcar = XDATCAR(arg)
    
    atoms = single(xdatcar, 1, line=[1], name="Center atom")
    atoms.put(single(xdatcar, 1, line=[2], show_info_flag=False, name="Collinear atom"))
    atoms.put(single(xdatcar, 2, line=[3], show_info_flag=False, name="Plane atoms"))
    
    xdatcar.read_all(atoms=atoms)
    dihe_angle = wrap(xdatcar, type="dh")
    
    line, factor = smooth(dihe_angle)
    atoms = atoms.get(True)
    cen_atom, collin_atom, plane_atom_1, plane_atom_2 = atoms[0], atoms[1], atoms[2], atoms[3]
    
    filename = (f"{cen_atom}_{collin_atom}_{plane_atom_1}_{plane_atom_2}" if arg.args.output == "" else arg.args.output)
    filename = arg.same_name(getcwd(), filename+".csv", factor)
    
    write_csv(filename, line, tuple(range(xdatcar.final_step)), f"{cen_atom}_{collin_atom}_{plane_atom_1}_{plane_atom_2}")
    print("Done!")

