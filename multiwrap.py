from os import getcwd
from numpy import array, newaxis
from copy import deepcopy
from kit.fundamental import AtomStep_Trj, Atom
from kit.vasp import XDATCAR
from kit.function import Distance
from kit.interface import couple, args, smooth
from kit.etc import write_csv

class Atom_Repeat(Atom):
    def __init__(self, input_obj=None, put_flag=False):
        Atom.__init__(self, input_obj, put_flag)
    def _in_list_check(self, atom):
        self._atom_list.append(atom)

if __name__ == "__main__":
    arg = args({"output": "distances"}, "XDATCAR")
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
    
    distance = Distance(xdatcar)
    cen_atom_step, mea_atom_step = AtomStep_Trj(xdatcar), AtomStep_Trj(xdatcar)
    center_positions, measure_positions = [], []
    center_atoms, measure_atoms = Atom_Repeat(xdatcar), Atom_Repeat(xdatcar)
    for atom_list, pos_list in zip(atom_lists, pos_lists):
        center_atoms.put(atom_list[0])
        measure_atoms.put(atom_list[1])
        center_positions.append(pos_list[0])
        measure_positions.append(pos_list[1])
    center_positions, measure_positions = array(center_positions).transpose((1, 0, 2)), array(measure_positions).transpose((1, 0, 2))
    cen_atom_step.steps, cen_atom_step.atoms, cen_atom_step.fractional_position = xdatcar.steps, center_atoms, center_positions
    mea_atom_step.steps, mea_atom_step.atoms, mea_atom_step.fractional_position = xdatcar.steps, measure_atoms, measure_positions
    distance.center_atom_step, distance.measure_atom_step = cen_atom_step, mea_atom_step
    distance.auto()
    distances_arr = distance.distance()
    lines, factor = smooth(distances_arr)
    arg.args.output = arg.same_name(getcwd(), arg.args.output, factor)
    columns = [f"{atom_list[0].get(True)[0]} {atom_list[1].get(True)[0]}" for atom_list in atom_lists]
    write_csv(arg.args.output, lines, xdatcar.steps.get(True), columns)
    print("Done!")
