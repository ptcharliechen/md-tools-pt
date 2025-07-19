from os import mkdir, path
from numpy import array, zeros
from kit.fundamental import Molecule
from kit.vasp import POSCAR
from kit.interface import args, lines

if __name__ == "__main__":
    arg = args({"output": "DOS", "fermi": [1, int, "fermi level correction. default is open (i.e. 1)"], "molecules": [None, str, "specify a file with molecule information"], "threshold": [0.0002, float]}, "vasprun.xml")
    
    poscar = POSCAR(arg)
    poscar.args.input = "POSCAR"
    atom_list = lines(poscar)
    
    orbital_order = ["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "dx2-y2"]
    while(1):
        orbitals = []
        tmp = input("Orbitals (s/px/py/pz/p/dxy/dxz/dyz/dz2/dx2-y2/d/all) [all]: ")
        if tmp == "" or tmp.lower() == "all":
            orbitals = [i for i in range(9)]
            break
        for i in tmp.split(","):
            if i.lower() == "p":
                orbitals.extend([1, 2, 3])
            elif i.lower() == "d":
                orbitals.extend([4, 5, 6, 7, 8])
            elif i.lower() not in orbital_order:
                print(f"Warning: '{i}' argument is an invalid argument.")
                break
            else:
                orbitals.append(orbital_order.index(i.lower()))
        else:
            break

    molecules = Molecule(arg)
    mkdir(arg.args.output)
    
    state, eigenvalue = [], []
    atoms_states, atoms_eigens = [], []
    with open(arg.args.input) as read_file:
        num_flag = False
        num = -2
        for line in read_file:
            if "efermi" in line:
                fermi_energy = float(line.split()[2]) if arg.args.fermi else 0.0
            elif "comment" in line and "ion" in line:
                num = -2
                num_flag = True
            elif num_flag and "</set>" in line:
                num_flag = False
                atoms_eigens.append(eigenvalue)
                atoms_states.append(state)
                state, eigenvalue = [], []
            elif num > -1 and num_flag:
                sp = line.split()
                eigenvalue.append(float(sp[1])-fermi_energy)
                state.append([float(i) for i in sp[2:11]])
            if num_flag:
                num += 1
    atoms_eigens = array(atoms_eigens)
    atoms_states = array(atoms_states)

    with open(path.join(arg.args.output, "list.csv"), "w") as write_file:
        write_file.write("No.,molecule,HOMO 1,HOMO 2,LUMO 1,LUMO 2,Gap,Step,Atoms\n")
    for mol in atom_list:
        atoms = mol.get()
        for key, molecule_list in molecules.molecule_dictionary.items():
            if atoms[0] in molecule_list:
                break
        DOS_orbitals = zeros((atoms_states.shape[1], len(orbitals)))
        DOS_total = zeros((atoms_states.shape[1], 9))
        for i in atoms:
            for idx, j in enumerate(orbitals):
                DOS_orbitals[:, idx] += atoms_states[i, :, j]
            DOS_total += atoms_states[i]
        with open(path.join(arg.args.output, molecules.molecule_kind[key]+".csv"), "w") as write_file:
            HOMO = zeros(2)
            LUMO = zeros(2)
            LUMO_flag = False
            for i in range(atoms_states.shape[1]):
                write_file.write(f"{atoms_eigens[key][i]:.4f},{sum(DOS_orbitals[i]):.4f}\n")
                if i > 0 and not LUMO_flag:
                    if sum(DOS_total[i]) < arg.args.threshold and sum(DOS_total[i-1]) > arg.args.threshold and atoms_eigens[key][i] < 0:
                        HOMO[0], HOMO[1] = atoms_eigens[key][i-1], atoms_eigens[key][i]
                    elif sum(DOS_total[i]) > arg.args.threshold and sum(DOS_total[i-1]) < arg.args.threshold and atoms_eigens[key][i] > 0:
                        LUMO[0], LUMO[1] = atoms_eigens[key][i-1], atoms_eigens[key][i]
                        LUMO_flag = True
        with open(path.join(arg.args.output, "list.csv"), "a") as list_file:
            list_file.write(f"{key+1},{molecules.molecule_kind[key]},{HOMO[0]:.4f},{HOMO[1]:.4f},{LUMO[0]:.4f},{LUMO[1]:.4f},{LUMO[1]-HOMO[0]:.4f},{LUMO[1]-LUMO[0]:.4f},{str(molecules.molecule_dictionary[key]).replace(',', ' ')}\n")
    print("Done!")
