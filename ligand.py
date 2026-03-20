from kit.args import args
arg = args({"output": "ligand", "mode": ["a", str, "atom mode (a), molecule mode (m), or cluster mode (c)"], \
            "bond": [None, str, "specify a file with user-defined bond length threshold"], \
            "dative": [None, float, "cutoff distance for dative bond"], \
            "scaling": [2.5, float, "set larger if you have large molecules"]}, "POSCAR", molecules=True, output_isdir=True)
from os import mkdir, path
from collections import Counter
from re import compile
from numpy import array, where
import numpy as np
from kit.fundamental import Atom, AtomStep_Single_Point
from kit.vasp import POSCAR
from kit.interface import blocks, read_bond
from kit.function import Ligand, HCE

if __name__ == "__main__":
    if arg.args.mode.lower() not in ["a", "m", "c"]:
        print("Warning: mode must be \"a\" (atom mode), \"m\" (molecule mode), or \"c\" (cluster mode). Reset as \"a\".")
        arg.args.mode = "a"
    elif not arg.args.mode.islower():
        arg.args.mode = arg.args.mode.lower()

    poscar = POSCAR(arg)
    
    atom_list, atom_list_flatten, atom_info = blocks(poscar, atom=Atom(), items=1, block_number=3 if poscar.args.mode == "c" else 2, print_output_flag=True, flatten_flag=True)
    atom_list, atom_info = atom_list[0], atom_info[0]

    poscar.read_all(atoms=atom_list_flatten)

    cen_atom, cluster_atom = Atom(poscar), Atom(poscar)
    if poscar.args.mode.lower() == "a" or poscar.args.mode.lower() == "m":
        cen_atom.put(atom_list[0])
        cluster_atom.put(atom_list[1])
    elif poscar.args.mode == "c":
        cen_atom.put(atom_list[0]); cen_atom.put(atom_list[1])
        cluster_atom.put(atom_list[0]); cluster_atom.put(atom_list[1]); cluster_atom.put(atom_list[2])

    center_atom_step, cluster_atom_step = AtomStep_Single_Point(poscar), AtomStep_Single_Point(poscar)
    center_pos, cluster_pos = [], []
    center_atom_step.atoms = cen_atom

    if poscar.args.mode == "c":
        for atom, pos in zip(atom_list_flatten.get(), poscar.fractional_position):
            if atom in cen_atom.get():
                center_pos.append(pos)
            cluster_pos.append(pos)
    else:
        for atom, pos in zip(atom_list_flatten.get(), poscar.fractional_position):
            if atom not in cen_atom.get():
                cluster_pos.append(pos)
            else:
                center_pos.append(pos)
    cluster_atom_step.atoms = cluster_atom
    center_atom_step.lock(); cluster_atom_step.lock()
    center_atom_step.fractional_position, cluster_atom_step.fractional_position = center_pos, cluster_pos
    
    cluster = HCE(poscar) if poscar.args.mode == "c" else Ligand(poscar)
    if cluster.args.dative is None or compile(r"^\d+\.?\d*$").match(str(cluster.args.dative)) is None:
        while(1):
            inp = input("Dative bond length [2.5]: ")
            if inp == "":
                cluster.args.dative = cluster.dative_bond = 2.5
                break
            elif compile(r"^\d+\.?\d*$").match(str(inp)) is not None:
                cluster.args.dative = cluster.dative_bond = float(inp)
                break
            else:
                print("Warning: Input a positive number.")
    else:
        cluster.dative_bond = cluster.args.dative
    
    cluster.center_atom_step, cluster.measure_atom_step = center_atom_step, cluster_atom_step
    
    if poscar.args.molecules is not None and path.isfile(poscar.args.molecules):
        cluster.read_molecules(poscar.args.molecules)
    elif poscar.args.molecules is not None:
        print(f"Warning: {poscar.args.molecules} does not exist.")
    
    read_bond(poscar)
    
    cluster.build_graph(molecule=True if cluster.molecules.molecule_dictionary != {} else False)
    cluster.graph_to_molecule_position()

    if poscar.args.mode == "c":
        cluster.run(atom_list[0], atom_list[1])
        anion_mole, cation_mole = [], []
        for mol, atoms in cluster.molecules.molecule_dictionary.items():
            for atom in atoms:
                if atom in atom_list[0].get():
                    cation_mole.append(mol)
                    break
                elif atom in atom_list[1].get():
                    anion_mole.append(mol)
                    break
    else:
        cluster.run()
    
    mkdir(poscar.args.output)
    mkdir(path.join(poscar.args.output, "clusters"))

    mole_kind = [kind.split('_')[0] for kind in cluster.molecules.molecule_kind]
    mole_kind_sorted = sorted(list(set(mole_kind)))

    idx_list = poscar.atom_step.atoms.index_list
    cluster_atoms_dict, cluster_position_dict = cluster.cluster_atoms, cluster.cluster_position
    for cen_mol, atoms in cluster_atoms_dict.items():
        if cluster.molecules.molecule_kind != []:
            mole_kind, cluster_type = cluster.cluster_types()
            if sum(cluster_type[cen_mol]) == 1:
                continue
            filename = "_".join([f"{kind}{type}" for (kind, type) in zip(mole_kind, cluster_type[cen_mol])])
        else:
            filename = cen_mol
        
        elements = []
        for atom in atoms:
            elements.append(poscar.elements[idx_list[atom]])
        
        with open(path.join(poscar.args.output, "clusters", f"{filename}.xyz"), "a") as write_file:
            write_file.write(f"{len(elements)}\n")
            write_file.write("center")
            for atom in cluster.molecules.molecule_dictionary[cen_mol]:
                write_file.write(f",{atom+poscar.args.atom}")
            write_file.write(",,bonded")

            for atom in cluster_atoms_dict[cen_mol][len(cluster.molecules.molecule_dictionary[cen_mol]):]:
                write_file.write(f",{atom+poscar.args.atom}")
            write_file.write("\n")
            for element, pos in zip(elements, cluster_position_dict[cen_mol]):
                write_file.write(f"{element}\t{pos[0]:.8f}\t{pos[1]:.8f}\t{pos[2]:.8f}\n")
    
    mole_pos = cluster.molecules.molecule_position
    with open(path.join(poscar.args.output, "ligand.csv"), "w") as write_file:
        write_file.write(f"center atom,element,molecule,bonded atom,element,molecule,distance")
        if poscar.args.mode == "c":
            write_file.write(",cluster type\n")
        else:
            write_file.write("\n")
        
        recorded_cen_moles = []
        if poscar.args.mode == "c":
            hce_cluster_type = {}
        neighbors = cluster.build_neighbors(cluster.cluster_graph)
        c2f, f2c = cluster.atom_step.c2f, cluster.atom_step.f2c
        
        for cen_mol in sorted(cluster.cluster_graph.keys()):
            if cen_mol in recorded_cen_moles:
                continue
            recorded_cen_moles.append(cen_mol)
            cen_atoms_idx = cluster.molecules.molecule_dictionary[cen_mol]
            cen_atoms_len = len(cen_atoms_idx)
            if cluster.args.mode != "c" and (cen_atoms_len == len(cluster_atoms_dict[cen_mol]) or cluster_atoms_dict[cen_mol] == []):
                continue
            elif cluster.args.mode == "c" and cen_mol not in cation_mole:
                continue
            elif cluster.args.mode == "c" and (cen_atoms_len == len(cluster_atoms_dict[cen_mol]) or cluster_atoms_dict[cen_mol] == []):
                cluster_atoms_dict[cen_mol].extend(cluster.molecules.molecule_dictionary[cen_mol])
                cluster_position_dict[cen_mol].extend(mole_pos[cen_mol])
                recorded_neighbor_moles = []
                for neighbor in neighbors[cen_mol]:
                    if neighbor in recorded_neighbor_moles:
                        continue
                    recorded_neighbor_moles.append(neighbor)
                    cluster_atoms_dict[cen_mol].extend(cluster.molecules.molecule_dictionary[neighbor])
                    mea_pos = mole_pos[neighbor]
                    dist_mat, period_image_mat = cluster.distance_matrix_cutoff(mole_pos[cen_mol], mea_pos, cluster.lattice, cluster.period_images, cluster.dative_bond, True)
                    for ii, jj in zip(*np.where(dist_mat < cluster.dative_bond)):
                        mea_pos[jj] += cluster.period_images[period_image_mat[ii, jj]]
                        mea_pos = cluster.image_shift(mea_pos, cluster.lattice, benchmark=jj, cartesian=0)
                        break
                    cluster_position_dict[cen_mol].extend(mea_pos)
                cluster_position_dict[cen_mol] = f2c(array(cluster_position_dict[cen_mol]), cluster.lattice)
            cen_pos = c2f(array(cluster_position_dict[cen_mol][:cen_atoms_len]), cluster.lattice)
            mea_pos = c2f(array(cluster_position_dict[cen_mol][cen_atoms_len:]), cluster.lattice)
            dist_vec = cluster.distance_matrix_func(cen_pos, mea_pos, cluster.lattice)
            
            for i, j in zip(*where(dist_vec < cluster.dative_bond)):
                cen_idx, mea_idx = cluster_atoms_dict[cen_mol][i], cluster_atoms_dict[cen_mol][cen_atoms_len+j]
                flag = False
                for mol, atoms in cluster.molecules.molecule_dictionary.items():
                    if cen_idx in atoms and mea_idx in atoms:
                        flag = True
                        break
                    if cen_idx in atoms:
                        cen_mol_num = mol
                    elif mea_idx in atoms:
                        mea_mol_num = mol
                if flag:
                    continue

                write_file.write(f"{cen_atoms_idx[i]+poscar.args.atom},{poscar.elements[idx_list[cen_atoms_idx[i]]]},")
                if cluster.molecules.molecule_kind != []:
                    write_file.write(f"{cluster.molecules.molecule_kind[cen_mol]},")
                else:
                    write_file.write(f"{cen_mol+1},")
                write_file.write(f"{cluster_atoms_dict[cen_mol][cen_atoms_len+j]+poscar.args.atom},{poscar.elements[idx_list[cluster_atoms_dict[cen_mol][cen_atoms_len+j]]]},")
                atom_idx = idx_list[cluster_atoms_dict[cen_mol][cen_atoms_len+j]]
                for ligand in cluster.ligand_molecules[cen_mol]:
                    if atom_idx in cluster.molecules.molecule_dictionary[ligand]:
                        if cluster.molecules.molecule_kind != []:
                            write_file.write(f"{cluster.molecules.molecule_kind[ligand]}")
                            break
                        else:
                            write_file.write(f"{ligand+1}")
                            break
                write_file.write(f",{dist_vec[i, j]:.4f}")
                if poscar.args.mode == "c":
                    anion_cation_counter = Counter()
                    recorded_ani_cat_mole = []
                    for mol in cluster.ligand_molecules[cen_mol]:
                        if mol in anion_mole and mol not in recorded_ani_cat_mole:
                            anion_cation_counter["anion"] += 1
                        elif mol in cation_mole and mol not in recorded_ani_cat_mole:
                            anion_cation_counter["cation"] += 1
                        recorded_ani_cat_mole.append(mol)
                    if anion_cation_counter["anion"] == 0:
                        hce_cluster_type[cen_mol] = "SSIP"
                        write_file.write(",SSIP")
                    elif anion_cation_counter["cation"] == 1:
                        hce_cluster_type[cen_mol] = "CIP"
                        write_file.write(",CIP")
                    elif anion_cation_counter["cation"] > 1:
                        hce_cluster_type[cen_mol] = "AGG"
                        write_file.write(",AGG")
                write_file.write("\n")
    
    if cluster.molecules.molecule_kind != []:
        with open(path.join(poscar.args.output, "cluster.csv"), "w") as write_file:
            cluster_type_count = Counter()
            if poscar.args.mode == "c":
                atoms_indices = f"atoms,{atom_info[0].replace(',', ';')} {atom_info[1].replace(',', ';')} {atom_info[2].replace(',', ';')}\n\n"
            else:
                atoms_indices = f"atoms,{atom_info[0].replace(',', ';')} {atom_info[1].replace(',', ';')}\n\n"
            if poscar.args.mode == "a":
                write_file.write("mode,atom\n")
            elif poscar.args.mode == "m":
                write_file.write("mode,molecule\n")
            elif poscar.args.mode == "c":
                write_file.write("mode,cluster\n")
            write_file.write(atoms_indices)
            for mol, ligand_count in cluster_type.items():
                if ligand_count is None:
                    continue
                elif arg.args.mode == "c" and mol in anion_mole:
                    continue
                ligand_count_str = ",".join([str(count) for count in ligand_count])
                cluster_type_count[ligand_count_str] += 1
            
            kinds = ",".join([str(_) for _ in mole_kind])
            write_file.write(f"{kinds},sum,percentage\n")
            for key, value in cluster_type_count.items():
                write_file.write(f"{key},{value},{value / sum(cluster_type_count.values()):.4f}\n")
            write_file.write('\n')

            if poscar.args.mode == "c":
                hce_cluster_counter = Counter(hce_cluster_type.values())
                atoms_sum = len(hce_cluster_type.values())
                write_file.write(f"SSIP,CIP,AGG\n")
                write_file.write(f"{hce_cluster_counter['SSIP']},{hce_cluster_counter['CIP']},{hce_cluster_counter['AGG']}\n")
                write_file.write(f"{hce_cluster_counter['SSIP']/atoms_sum:.4f},{hce_cluster_counter['CIP']/atoms_sum:.4f},{hce_cluster_counter['AGG']/atoms_sum:.4f}\n\n")

            write_file.write(f"center molecule,{kinds}\n")
            for key in sorted(cluster_type.keys()):
                if poscar.args.mode == "c" and key not in cation_mole:
                    continue
                write_file.write(f"{cluster.molecules.molecule_kind[key]}")
                for i in cluster_type[key]:
                    write_file.write(f",{i}")
                write_file.write("\n")

    if poscar.args.mode == "c":
        poscar_first_solvation = POSCAR(arg)
        poscar_first_solvation.args.output = path.join(poscar.args.output, "POSCAR_first_shield")
        poscar_first_solvation.title = "POSCAR_first_shield"
        recorded_moles, positions, elements = [], [], []
        for moles in cluster.ligand_molecules.values():
            for mole in moles:
                if mole not in recorded_moles:
                    recorded_moles.append(mole)
                    for atom in cluster.molecules.molecule_dictionary[mole]:
                        elements.append(poscar.elements[idx_list[atom]])
                    positions.extend(mole_pos[mole])
        poscar_first_solvation.lattice = poscar.lattice
        poscar_first_solvation.rearrange(elements, positions)
        poscar_first_solvation.write_all()
    print("Done!")
