from kit.args import args
arg = args({"output": "distances"}, "XDATCAR")
from os import getcwd
from copy import deepcopy
import numpy as np
from kit.fundamental import AtomStep_Trj, Atom
from kit.vasp import XDATCAR
from kit.function import Distance
from kit.interface import couple, smooth
from kit.etc import write_csv

class Atom_Repeat(Atom):
    def __init__(self, input_obj=None, put_flag=False):
        Atom.__init__(self, input_obj, put_flag)
    def _in_list_check(self, atom):
        self._atom_list.append(atom)

if __name__ == "__main__":
    xdatcar = XDATCAR(arg)
    
    atom_lists, total_lists = couple(xdatcar, 2, atom=Atom_Repeat(), Atom_split_flag=True, flatten_flag=True)

    xdatcar.read_all(atoms=total_lists)
    pos_lists = []
    frac_pos = xdatcar.fractional_position
    total_lists = total_lists.get()
    for atom_list in atom_lists:
        pos = []
        for atom in atom_list:
            pos.append(frac_pos[:, total_lists.index(atom.get()[0])])
        pos_lists.append(deepcopy(pos))
    
    distances = []
    for atom_list, pos_list in zip(atom_lists, pos_lists):
        distance = Distance(xdatcar)
        cen_atom_step, mea_atom_step = AtomStep_Trj(xdatcar), AtomStep_Trj(xdatcar)
        cen_atom_step.steps, cen_atom_step.atoms, cen_atom_step.fractional_position = xdatcar.steps, atom_list[0], pos_list[0][np.newaxis, :].transpose((1, 0, 2))
        mea_atom_step.steps, mea_atom_step.atoms, mea_atom_step.fractional_position = xdatcar.steps, atom_list[1], pos_list[1][np.newaxis, :].transpose((1, 0, 2))
        distance.center_atom_step, distance.measure_atom_step = cen_atom_step, mea_atom_step
        distance.auto()
        distances.append(distance.distance())
    
    lines, factor = smooth(np.array(distances)[:, 0, :])
    arg.args.output = arg.same_name(getcwd(), arg.args.output, factor)
    columns = [f"{atom_list[0].get(True)[0]} {atom_list[1].get(True)[0]}" for atom_list in atom_lists]
    write_csv(arg.args.output, lines, xdatcar.steps.get(True), columns)
    print("Done!")
