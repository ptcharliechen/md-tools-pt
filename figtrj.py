from kit.args import args
arg = args({"output": "POSCAR"}, "XDATCAR")
from kit.fundamental import Atom, Step
from kit.vasp import XDATCAR, POSCAR
from kit.interface import line, step

if __name__ == "__main__":
    xdatcar = XDATCAR(arg)
    
    atom_list = line(xdatcar)
    step_list = step(xdatcar)
    
    elements_complement = []
    atom_list_complement = Atom(xdatcar)
    for atom in range(len(xdatcar.elements)):
        if atom not in atom_list.get():
            atom_list_complement.put(atom)
            elements_complement.append(xdatcar.elements[atom])
    while(1):
        tmp = input("\nImport background atoms (y/n) [y]: ")
        if tmp == "":
            tmp = "y"
        if tmp.lower() == "y":
            tmp = input("Background step [1]: ")
            if tmp == "":
                tmp = "1"
            if tmp.isdigit():
                background = Step(xdatcar)
                background.put(tmp)
                break
            else:
                print("Warning: Input error.")
        elif tmp.lower() == "n":
            background = -1
            break
        else:
            print("Warning: Input error.")
    xdatcar.read_all(step_list, atom_list)
    if background != -1:
        xdatcar_background = XDATCAR(arg)
        xdatcar_background.read_all(background, atom_list_complement)
    poscar = POSCAR(arg)
    poscar.args, poscar.title, poscar.elements, poscar.lattice = xdatcar.args, xdatcar.title, xdatcar.elements, xdatcar.lattice
    elements = xdatcar.elements
    ele_dict = {}
    frac_pos = xdatcar.fractional_position.transpose((1, 0, 2))
    for idx, atom in enumerate(atom_list.get()):
        if elements[atom] not in ele_dict.keys():
            ele_dict[elements[atom]] = list(frac_pos[idx])
        else:
            ele_dict[elements[atom]] += list(frac_pos[idx])
    if background != -1:
        frac_pos_background = xdatcar_background.fractional_position
        for idx, atom in enumerate(atom_list_complement.get()):
            if elements[atom] not in ele_dict.keys():
                ele_dict[elements[atom]] = [frac_pos_background[0][idx]]
            else:
                ele_dict[elements[atom]].append(frac_pos_background[0][idx])
    poscar_elements, positions = [], []
    for element, position in ele_dict.items():
        for _ in range(len(position)):
            poscar_elements.append(element)
        positions += position
    poscar.atoms.put(range(len(poscar_elements)))
    poscar.atoms.elements = poscar_elements
    poscar.fractional_position = positions
    poscar.write_all()
    print("Done!")
