from kit.args import args
arg = args({"output": "reaction", "cpu": [3, int, "specify the maximum number of applied CPU cores"], \
            "samename": [1, int, "0: without same name check, 1: with same name check"], "bond": [None, str, "specify a file with bond threshold"], \
            "scaling": [2.5, float, "set larger if you have large molecules"]}, "XDATCAR", molecules=True, output_isdir=True)
from os import mkdir, path, cpu_count
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from kit.fundamental import Atom, Step, Graph, AtomStep_Trj
from kit.interface import step, read_bond
from kit.function import Reaction
from kit.vasp import XDATCAR

def function(num, atom_step):
    print(f"Calculating Part {num}...")
    
    reaction_init = Graph()
    reaction_init.args = atom_step.args
    
    reaction_init.atom_step = atom_step
    reaction_init.build_graph(wrap=False if num == "1" else True, molecule=True if num == "1" else False)
    reaction_init.graph_to_molecule_position(reaction_init.adjacent_matrix)

    reaction = Reaction(reaction_init)
    reaction.period_images = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0],
                                       [0, 1, 0], [0, -1, 0], [1, 1, 0],
                                       [1, -1, 0], [-1, 1, 0], [-1, -1, 0]])
    reaction.adjacent_matrix = reaction_init.adjacent_matrix

    reaction.run()
    
    return reaction.bond_info

if __name__ == "__main__":
    xdatcar = XDATCAR(arg)
    
    step_list = step(xdatcar)
    read_bond(xdatcar)
    
    xdatcar.read_molecules(xdatcar.args.molecules)
    atom_list = Atom(xdatcar)
    molecules = xdatcar.atom_step.molecules
    for mole_list in molecules.molecule_dictionary.values():
        for atom in mole_list:
            atom_list.put(atom)
    atom_list.sort()
    
    print("Reading XDATCAR...")
    xdatcar.read_all(steps=step_list, atoms=atom_list)

    if arg.args.cpu <= 0:
        arg.args.cpu = 3 if cpu_count() > 4 else 1
    elif arg.args.cpu > cpu_count():
        arg.args.cpu = cpu_count()

    positions = np.array_split(xdatcar.atom_step.fractional_position, arg.args.cpu)
    if step_list == "all":
        step_list = xdatcar.atom_step.steps
    steps = np.array_split(step_list.get(), arg.args.cpu)
    atom_steps = []
    nums = []
    for idx, (_step, position) in enumerate(zip(steps, positions)):
        atom_step = AtomStep_Trj(xdatcar)
        atom_step.molecules.molecule_kind = []
        atom_step.molecules.molecule_dictionary = defaultdict(list)
        atom_step.molecules.molecule_position = defaultdict(list)
        atom_step.atoms.elements = [xdatcar.elements[element] for element in atom_list.get()]
        if idx > 0:
            _step_list = [steps[idx-1][-1]] + list(_step)
            _step = Step(xdatcar)
            atom_step.args.step = _step.args.step = max(steps[idx-1])
            _step.put(_step_list)
            atom_step.steps = _step
            atom_step.fractional_flag = True
            position = list(position)
            position.insert(0, positions[idx-1][-1])
            atom_step.fractional_position = np.array(position)
        elif idx == 0:
            _step_list = list(_step)
            _step = Step(xdatcar)
            _step.put(_step_list)
            atom_step.steps = _step
            atom_step.fractional_flag = True
            atom_step.fractional_position = position
        nums.append(str(idx+1))
        atom_steps.append(atom_step)
    
    with ProcessPoolExecutor(max_workers=arg.args.cpu) as executor:
        futures = [executor.submit(function, num, atom_step) for num, atom_step in zip(nums, atom_steps)]
        results = [future.result() for future in futures]
    
    bond_info = {1: [], 2: [], 3: []}
    for result in results:
        bond_info[1].extend(result[1])
        bond_info[2].extend(result[2])
        bond_info[3].extend(result[3])
    
    elements = xdatcar.atom_step.atoms.elements
    idx_list = xdatcar.atom_step.atoms.index_list
    if path.isdir(arg.args.output) and arg.args.samename:
        rmdir(arg.args.ouput)
    mkdir(arg.args.output)
    with open(path.join(arg.args.output, "bond_breaking.dat"), "w") as write_file:
        if bond_info[1] + bond_info[2] == []:
            write_file.write("No bond is broken\n")
        else:
            for i, mole_list in molecules.molecule_dictionary.items():
                counter = 0
                write_file.write(f"Molecule: {molecules.molecule_kind[i]}\n")
                removes = []
                for j, info in enumerate(bond_info[1]):
                    if info[1] in mole_list or info[2] in mole_list:
                        if any([bond_info[1] == info[1] and bond_info[2] == info[2] for bond_info in bond_info[1][j+1:]]):
                            removes.append(info)
                            continue
                        write_file.write(f"\tStep: {info[0]+xdatcar.args.step}, New fragments: No\n")
                        write_file.write(f"\t\tAtom {info[1]}, Element: {elements[idx_list[info[1]]]}, Position: ({info[3][0]:.2f}, {info[3][1]:.2f}, {info[3][2]:.2f})\n")
                        write_file.write(f"\t\tAtom {info[2]}, Element: {elements[idx_list[info[2]]]}, Position: ({info[4][0]:.2f}, {info[4][1]:.2f}, {info[4][2]:.2f})\n")
                        write_file.write(f"\t\t{elements[idx_list[info[1]]]}-{elements[idx_list[info[2]]]} threshold: {molecules.threshold[f'{elements[idx_list[info[1]]]}-{elements[idx_list[info[2]]]}']:.2f} A. Distance: {info[5]:.2f} A\n")
                        counter += 1
                        removes.append(info)
                if removes != []:
                    for remove in removes:
                        bond_info[1].remove(remove)
                removes = []
                for j, info in enumerate(bond_info[2]):
                    if info[1] in mole_list or info[2] in mole_list:
                        for build_info in bond_info[3]:
                            if info[1] == build_info[1] and info[2] == build_info[2] and build_info[0] > info[0]:
                                bond_info[3].remove(build_info)
                                break
                        else:
                            if any([bond_info[1] == info[1] and bond_info[2] == info[2] for bond_info in bond_info[2][j+1:]]):
                                removes.append(info)
                                continue
                            write_file.write(f"\tStep: {info[0]+xdatcar.args.step}, New fragments: Yes\n")
                            write_file.write(f"\t\tAtom {info[1]}, Element: {elements[idx_list[info[1]]]}, Position: ({info[3][0]:.2f}, {info[3][1]:.2f}, {info[3][2]:.2f}), Fragment: {info[6]}\n")
                            write_file.write(f"\t\tAtom {info[2]}, Element: {elements[idx_list[info[2]]]}, Position: ({info[4][0]:.2f}, {info[4][1]:.2f}, {info[4][2]:.2f}), Fragment: {info[7]}\n")
                            write_file.write(f"\t\t{elements[idx_list[info[1]]]}-{elements[idx_list[info[2]]]} threshold: {molecules.threshold[f'{elements[idx_list[info[1]]]}-{elements[idx_list[info[2]]]}']:.2f} A. Distance: {info[5]:.2f} A\n")
                            counter += 1
                        removes.append(info)
                if removes != []:
                    for remove in removes:
                        bond_info[2].remove(remove)
                write_file.write(f"Break {counter} bond")
                if counter > 1:
                    write_file.write("s")
                write_file.write("\n\n")
    with open(path.join(arg.args.output, "bond_building.dat"), "w") as write_file:
        if bond_info[3] == []:
            write_file.write("No bond is built\n")
        else:
            for i, mole_list in molecules.molecule_dictionary.items():
                counter = 0
                write_file.write(f"Molecule: {molecules.molecule_kind[i]}\n")
                removes = []
                for j, info in enumerate(bond_info[3]):
                    if info[1] in mole_list or info[2] in mole_list:
                        if any([bond_info[1] == info[1] and bond_info[2] == info[2] for bond_info in bond_info[3][j+1:]]):
                            removes.append(info)
                            continue
                        for j, inner_mole_list in molecules.molecule_dictionary.items():
                            if info[1] in inner_mole_list:
                                mole_min = molecules.molecule_kind[j]
                            if info[2] in inner_mole_list:
                                mole_max = molecules.molecule_kind[j]
                        write_file.write(f"\tStep: {info[0]+xdatcar.args.step}, New fragment: ")
                        if info[6] != "":
                            write_file.write(f"{info[6]}\n")
                        else:
                            write_file.write("No\n")
                        write_file.write(f"\t\tAtom {info[1]}, Element: {elements[idx_list[info[1]]]}, Position: ({info[3][0]:.2f}, {info[3][1]:.2f}, {info[3][2]:.2f}), Molecule: {mole_min}\n")
                        write_file.write(f"\t\tAtom {info[2]}, Element: {elements[idx_list[info[2]]]}, Position: ({info[4][0]:.2f}, {info[4][1]:.2f}, {info[4][2]:.2f}), Molecule: {mole_max}\n")
                        write_file.write(f"\t\t{elements[idx_list[info[1]]]}-{elements[idx_list[info[2]]]} threshold: {molecules.threshold[f'{elements[idx_list[info[1]]]}-{elements[idx_list[info[2]]]}']:.2f} A. Distance: {info[5]:.2f} A\n")
                        counter += 1
                        removes.append(info)
                if removes != []:
                    for remove in removes:
                        bond_info[3].remove(remove)
                write_file.write(f"Build {counter} bond")
                if counter > 1:
                    write_file.write("s")
                write_file.write("\n\n")
    print("Done!")
