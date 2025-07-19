from kit.fundamental import Atom, Step
from kit.vasp import XDATCAR, POSCAR
from kit.interface import args, line, step

if __name__ == "__main__":
    arg = args({"output": "POSCAR"}, "XDATCAR")
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
        tmp = input("\nImport Background atoms (y/n) [y]: ")
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
    xdatcar.read_step_atom(step_list, atom_list)
    if background != -1:
        xdatcar_2 = XDATCAR(arg)
        xdatcar_2.read_step_atom(background, atom_list_complement)
    poscar = POSCAR(arg)
    poscar.args, poscar.title, poscar.elements, poscar.lattice = xdatcar.args, xdatcar.title, xdatcar.elements, xdatcar.lattice
    elements = xdatcar.elements
    elemDict = {}
    fracPos = xdatcar.fractional_position.transpose((1, 0, 2))
    for idx, atom in enumerate(atom_list.get()):
        if elements[atom] not in elemDict.keys():
            elemDict[elements[atom]] = list(fracPos[:, idx])
        else:
            elemDict[elements[atom]] += list(fracPos[:, idx])
    if background != -1:
        fracPos_2 = xdatcar_2.fractional_position
        for idx, atom in enumerate(atom_list_complement.get()):
            if elements[atom] not in elemDict.keys():
                elemDict[elements[atom]] = [fracPos_2[idx][0]]
            else:
                elemDict[elements[atom]].append(fracPos_2[idx][0])
    poscar_elements, positions = [], []
    for element, position in elemDict.items():
        for _ in range(len(position)):
            poscar_elements.append(element)
        positions += position
    poscar.elements = poscar_elements
    poscar.fractional_position = positions
    poscar.write_all()
    print("Done!")
