from numpy import where
from collections import defaultdict
from kit.vasp import POSCAR
from kit.interface import line, args
from kit.function import Ligand
from molecule import threshold

if __name__ == "__main__":
    #Preprocess
    arg = args({"output": ["neighbors", "output file name includes atom numbers of every molecule"], \
                "bond": [None, str, "specify a file with user-defined bond length threshold"], \
                "molecules": [None, str, "molecule information file"]}, "POSCAR")
    #Get atoms
    poscar = POSCAR(arg)
    
    atoms = line(poscar)
    poscar.read_all(atoms=atoms)
    molecule = Ligand()
    molecule.bridge = poscar
    #Set the bond threshold
    if arg.args.bond is not None:
        molecule.read_bond()
    else:
        bondType = threshold()
        if bondType is not None:
            molecule.threshold = bondType
    molecule.build_molecule_position()
    molecule.build_graph()
    graph = molecule.adjacent_matrix
    poscar.read_molecules()
    MoleKinds, moleDict = poscar.atom_step.molecules.molecule_kind, poscar.atom_step.molecules.molecule_dictionary
    AtomKind, MoleculeElements = {}, defaultdict(list)
    ElementsList = list(set(molecule.elements))
    for num in sorted(moleDict.keys()):
        for atom in moleDict[num]:
            surrondingElements = defaultdict(int)
            for _ in where(graph[atom])[0]:
                if _ == atom:
                    continue
                surrondingElements[molecule.elements[_]] += 1
                element = molecule.elements[_]
            IdxList = poscar.atom_step.index_list
            AtomKind[atom] = [MoleKinds[num], molecule.elements[IdxList[atom]], str(dict(surrondingElements))]
            if molecule.elements[IdxList[atom]] not in MoleculeElements[MoleKinds[num]]:
                MoleculeElements[MoleKinds[num]].append(molecule.elements[IdxList[atom]])
    
    with open(arg.args.output+".csv", "w") as write_file:
        recorded = []
        for MoleKind in set(MoleKinds):
            write_file.write(f"Molecule: {MoleKind}\n")
            for element in MoleculeElements[MoleKind]:
                write_file.write(f"Element: {element}\n")
                for AtomI, ValI in AtomKind.items():
                    if AtomI in recorded:
                        continue
                    if ValI[0] == MoleKind and ValI[1] == element:
                        recorded.append(AtomI)
                        if ValI[2] == "{}":
                            write_file.write(f"None:\n{AtomI}\n")
                        else:
                            atoms = []
                            write_file.write(f"{ValI[2]}:\n")
                            atoms.append(int(AtomI))
                            for AtomJ, ValJ in AtomKind.items():
                                if AtomJ in recorded:
                                    continue
                                if ValJ[0] == MoleKind and ValJ[1] == element and ValJ[2] == ValI[2]:
                                    atoms.append(int(AtomJ))
                                    recorded.append(AtomJ)
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
