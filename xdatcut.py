from kit.args import args
arg = args({"output": "XDATCAR_1"}, "XDATCAR")
from kit.vasp import XDATCAR
from kit.interface import step_lines

if __name__ == "__main__":
    xdatcar = XDATCAR(arg)
    steps, steps_flatten = step_lines(xdatcar, flatten_flag=True)
    xdatcar.read_all(steps=steps_flatten)
    for idx, step in enumerate(steps):
        args_new = args({"output": f"XDATCAR_{idx+1}"}, "XDATCAR")
        xdatcar_new = XDATCAR(xdatcar)
        xdatcar_new.args = args_new
        xdatcar_new.atom_step.split(xdatcar.atom_step, step, xdatcar.atoms)
        xdatcar_new.write_all()
    print("Done!")

