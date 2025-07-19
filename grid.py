from re import compile
import numpy as np
from scipy.ndimage import zoom
from kit.interface import args, step, line
from kit.vasp import XDATCAR
from kit.accelerate import shift_to_origin

if __name__ == "__main__":
    arg = args({"output": "grid", "molecules": [None, str, "specify a file with molecule information"]}, "XDATCAR")
    xdatcar = XDATCAR(arg)
    atom_list = line(xdatcar)
    step_list = step(xdatcar)
    
    while(1):
        tmp = input("mesh size [10]: ")
        if tmp == "":
            mesh = 10
            break
        elif tmp.isdigit():
            mesh = int(tmp)
            break
        else:
            print("Warning: Input error.")
    
    while(1):
        tmp = input("projection direction [xy]: ")
        if tmp == "":
            direction = [0, 1]
            break
        elif len(tmp) != 2:
            print("Warning: Input error.")
            continue
        direction = []
        for c in tmp.lower():
            if c == 'x':
                direction.append(0)
            elif c == 'y':
                direction.append(1)
            elif c == 'z':
                direction.append(2)
            else:
                print("Warning: Input error.")
                break
        else:
            if 2 not in direction:
                direction.append(2)
                projected_dir = 'z'
            elif 1 not in direction:
                direction.append(1)
                projected_dir = 'y'
            elif 0 not in direction:
                direction.append(0)
                projected_dir = 'x'
            break
    
    while(1):
        tmp = input(f"lower bound of {projected_dir} direction [0]: ")
        if tmp == "":
            lower = 0
            break
        elif compile(r"\d+\.?\d*").match(tmp) is not None:
            lower = float(tmp)
            break
        else:
            print("Warning: Input error.")
    while(1):
        tmp = input(f"upper bound of {projected_dir} direction [1]: ")
        if tmp == "":
            upper = 1
            break
        elif compile(r"\d+\.?\d*").match(tmp) is not None:
            upper = float(tmp)
            break
        else:
            print("Warning: Input error.")
    
    xdatcar.read_all(atoms=atom_list, steps=step_list)
    positions_idx = np.floor(shift_to_origin(xdatcar.fractional_position) * mesh).astype(int)
    atom_density = np.zeros((mesh, mesh), dtype=float)
    for positions_step_idx in positions_idx:
        for position_idx in positions_step_idx:
            if lower < position_idx[direction[2]] < upper:
                atom_density[position_idx[direction[0]], position_idx[direction[1]]] += 1
    atom_density /= positions_idx.shape[0]
    tiled_atom_density = np.tile(atom_density, (3, 3))

    while(1):
        yn = input("\nWould you like to smooth the curve (y/n): ")
        if yn.lower() == 'y':
            tiled_atom_density = zoom(tiled_atom_density, zoom=10, order=3)
            tiled_atom_density[tiled_atom_density < 0] = 0
            start, end = 10 * mesh, 20 * mesh + 1
            break
        elif yn.lower() == 'n':
            start, end = mesh, 2 * mesh + 1
            break
        else:
            print("Warning: Input 'y' or 'n'.")
    center_atom_density = tiled_atom_density[start:end, start:end]

    np.savetxt(arg.args.output+".csv", center_atom_density, delimiter=",", fmt="%.4f")

    print("Done!")
