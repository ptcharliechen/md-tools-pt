from kit.vasp import POSCAR
from kit.interface import line, args, threshold
from kit.fundamental import Graph

if __name__ == "__main__":
    #Preprocess
    arg = args({"output": ["molecules", "output file name includes atom numbers of every molecule"], \
                "bond": [None, str, "specify a file with user-defined bond length threshold"], \
                "wrap": [1, int, "whether to wrap the molecules"]}, "POSCAR")
    #Get atoms
    poscar = POSCAR(arg)
    
    atoms = line(poscar)
    poscar.read_all(atoms=atoms)
    molecules = Graph(poscar)
    #Set the bond threshold
    if arg.args.bond is not None:
        molecules.molecules.read_bond()
    else:
        bond_type = threshold()
        if bond_type is not None:
            molecules.threshold = bond_type
    molecules.build_graph(wrap=arg.args.wrap, molecule=False)
    molecules.graph_to_molecule_position(wrap=arg.args.wrap)
    #Write the file
    with open(f"{arg.args.output}.csv", "w") as write_file:
        for mole in molecules.molecules.molecule_dictionary.values():
            for idx, atom in enumerate(mole):
                write_file.write(f"{atom}" if idx == 0 else f",{atom}")
            write_file.write("\n")
    print("Done!")
