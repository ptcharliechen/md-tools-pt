from kit.args import args
arg = args({"output": "direction"}, "XDATCAR")
from kit.vasp import XDATCAR
from kit.accelerate import image_shift
from kit.interface import line, smooth
from kit.etc import write_csv

if __name__ == "__main__":
    xdatcar = XDATCAR(arg)
    
    atom_list = line(xdatcar)

    while(1):
        cart_dir = input("Which direction you would like to discover (x/y/z) [z]: ")
        if cart_dir.lower() == 'x' or cart_dir.lower() == 'y' or cart_dir.lower() == 'z':
            break
        elif cart_dir == "":
            cart_dir = 'z'
            break
        else:
            print("Warning: Input 'x', 'y', or 'z'.")

    xdatcar.read_all(atoms=atom_list)

    direction = image_shift(xdatcar.fractional_position, xdatcar.lattice)
    if cart_dir.lower() == 'x':
        pos = direction[:, :, 0]
    elif cart_dir.lower() == 'y':
        pos = direction[:, :, 1]
    elif cart_dir.lower() == 'z':
        pos = direction[:, :, 2]
    
    lines, factor = smooth(pos.T)
    if factor:
        arg.args.output += f"_Alpha={factor}"
    write_csv(arg.args.output, lines, list(range(pos.shape[0])), atom_list.get(True))
    print("Done!")

