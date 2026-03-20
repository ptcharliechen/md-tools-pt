from kit.args import args
#Preprocess
arg = args({"output": ["neighbors", "output file name includes atom numbers of every molecule"], \
            "bond": [None, str, "specify a file with user-defined bond length threshold"], \
            "mode": ["t", str, "type mode (t), or molecule mode (m)"]}, "POSCAR", molecules=True)
from numpy import where
from collections import defaultdict
from re import match
from kit.vasp import POSCAR
from kit.interface import line, read_bond
from kit.fundamental import Graph

if __name__ == "__main__":
    #Get atoms
    poscar = POSCAR(arg)
    
    atoms = line(poscar)
    poscar.read_all(atoms=atoms)
    molecule = Graph(poscar)

    #Set the bond threshold
    read_bond(poscar)
    poscar.read_molecules()
    molecule.build_graph()
    graph = molecule.adjacent_matrix
    
    if arg.args.mode.lower() == "m":
        mole_kinds, mole_dict = poscar.atom_step.molecules.molecule_kind, poscar.atom_step.molecules.molecule_dictionary
    else:
        mole_kinds, mole_dict = [_.split("_")[0] for _ in poscar.atom_step.molecules.molecule_kind], poscar.atom_step.molecules.molecule_dictionary
    atom_kind, molecule_elements = {}, defaultdict(list)
    
    ElementsList = list(set(molecule.elements))
    for num in sorted(mole_dict.keys()):
        for atom in mole_dict[num]:
            surrondingElements = defaultdict(int)
            for _ in where(graph[molecule.atoms.index_list[atom]])[0]:
                if _ == atom:
                    continue
                surrondingElements[molecule.elements[_]] += 1
                element = molecule.elements[_]
            IdxList = poscar.atom_step.atoms.index_list
            atom_kind[atom] = [mole_kinds[num], molecule.elements[IdxList[atom]], str(dict(surrondingElements))]
            if molecule.elements[IdxList[atom]] not in molecule_elements[mole_kinds[num]]:
                molecule_elements[mole_kinds[num]].append(molecule.elements[IdxList[atom]])

    with open(arg.args.output+".csv", "w") as write_file:
        recorded = []
        mole_kinds_sorted = sorted(set(mole_kinds)) if arg.args.mode.lower() == "t" else sorted(set(mole_kinds), key=lambda x: (match(r'(\D*)(\d+)', x).groups()[0], int(match(r'(\D*)(\d+)', x).groups()[1])))
        for mole_kind in mole_kinds_sorted:
            write_file.write(f"Molecule: {mole_kind}\n")
            for element in molecule_elements[mole_kind]:
                write_file.write(f"Element: {element}\n")
                for atom_i, val_i in atom_kind.items():
                    if atom_i in recorded:
                        continue
                    if val_i[0] == mole_kind and val_i[1] == element:
                        recorded.append(atom_i)
                        if val_i[2] == "{}":
                            write_file.write(f"None,total: 1 atoms\n{atom_i}\n")
                        else:
                            atoms = []
                            write_file.write(val_i[2].replace('{','').replace('}','').replace("'","")+",")
                            atoms.append(int(atom_i))
                            for atom_j, val_j in atom_kind.items():
                                if atom_j in recorded:
                                    continue
                                if val_j[0] == mole_kind and val_j[1] == element and val_j[2] == val_i[2]:
                                    atoms.append(int(atom_j))
                                    recorded.append(atom_j)
                            write_file.write(f"total: {len(atoms)} atoms\n")
                            for idx, i in enumerate(sorted(atoms)):
                                if idx == 0:
                                    write_file.write(f"{i}")
                                else:
                                    write_file.write(f",{i}")
                            write_file.write("\n")
                    else:
                        continue
            write_file.write("\n")
    print("Done!")

