from kit.vasp import XDATCAR
from kit.function import Mobility
from kit.interface import args, line, smooth
from kit.etc import write_csv

if __name__ == "__main__":
    arg = args({"output": "Mobility"}, "XDATCAR")

    while(1):
        tmp = input("Origin step: ")
        if not tmp.isdigit():
            print("Warning: Input a positive integer.")
        else:
            basis = int(tmp)
            break

    xdatcar = XDATCAR(arg)
    atom_list = line(xdatcar)
    print("Wait a moment!")
    xdatcar.read_all(atoms=atom_list)

    mobility = Mobility()
    mobility.bridge = xdatcar
    mobility.benchmark(basis)
    dist = mobility.mobility()

    lines, factor = smooth(dist)
    arg.args.output = arg.same_name(arg.args.output, factor)
    write_csv(arg.args.output, lines, list(range(xdatcar.final_step)), atom_list.get())
    print("Done!")
