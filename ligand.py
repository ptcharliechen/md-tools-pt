from os import mkdir, path
from re import compile
from numpy import array, where
import numpy as np
from kit.fundamental import Atom, AtomStep_Single_Point
from kit.vasp import POSCAR
from kit.interface import line, args, threshold
from kit.function import Ligand

class HCE(Ligand):
    def __init__(self, input_obj=None):
        Ligand.__init__(self, input_obj)
    def _build_cluster_graph(self, molecule_position, center_avg_molecules, measure_avg_molecules):
        self._cluster_adj_list = self.defaultdict(dict)
        threshold = self._lattice_bond * self._search_scaling
        lattice = self._AtomStep.lattice
        cluster_candidates = self.defaultdict(list)
        for mol, cen_avg_mole in zip(self._cen_mole, center_avg_molecules):
            if mol in self._cation_moles:
                mea_mole = self._mea_mole
                mea_pos = measure_avg_molecules
            else:
                mea_mole, mea_pos = [], []
                for idx, pos in zip(self._mea_mole, measure_avg_molecules):
                    if idx in self._cation_moles:
                        mea_mole.append(idx)
                        mea_pos.append(pos)
            
            dist_mat = self.distance_matrix_wrap(cen_avg_mole[np.newaxis, :], array(mea_pos), lattice, self._period_images)
            for _, j, _ in zip(*np.where(dist_mat < threshold)):
                cluster_candidates[mol].append(mea_mole[j])
        
        for cen_mole, mea_moles in cluster_candidates.items():
            cen_pos = molecule_position[cen_mole]
            for mol_id in mea_moles:
                if mol_id == cen_mole:
                    continue
                mea_pos = molecule_position[mol_id]
                dist_mat = self.distance_matrix_wrap(cen_pos, mea_pos, lattice, self._period_images)
                s = np.sum(dist_mat < self._lattice_bond)
                if s:
                    self._cluster_adj_list[mol_id][cen_mole] = self._cluster_adj_list[cen_mole][mol_id] = s
    def _graph_to_cluster_position(self, cluster_graph, molecule_dictionary, molecule_position):
        cluster_pos, cluster_atoms, ligand_molecules = self.defaultdict(list), self.defaultdict(list), self.defaultdict(list)
        
        recorded_cations = []
        for anion_mole in self._anion_moles:
            cluster_moles = self.walking_list(anion_mole, cluster_graph)
            cluster_avg_pos = []
            for cluster_mole in cluster_moles:
                cluster_atoms[anion_mole].extend(molecule_dictionary[cluster_mole])
                ligand_molecules[anion_mole].append(cluster_mole)
                if cluster_mole in cluster_graph[anion_mole].keys() and cluster_graph[anion_mole][cluster_mole] > 1:
                    for _ in range(cluster_graph[anion_mole][cluster_mole] - 1):
                        ligand_molecules[anion_mole].append(cluster_mole)
                mea_pos, mea_avg_pos = molecule_position[cluster_mole], np.average(molecule_position[cluster_mole], axis=0)
                if cluster_avg_pos != []:
                    dist_mat = self.distance_matrix_wrap(array(cluster_avg_pos), mea_avg_pos[np.newaxis, :], self._AtomStep.lattice, self._period_images)
                    _, _, k = np.unravel_index(np.argmin(dist_mat), dist_mat.shape)
                    mea_pos[0] += self._period_images[k]
                    mea_pos = self.image_shift(mea_pos, self._AtomStep.lattice, cartesian=0)
                cluster_pos[anion_mole].extend(list(mea_pos))
                cluster_avg_pos.append(mea_avg_pos)
                if cluster_mole != anion_mole and cluster_mole in cluster_atoms.keys():
                    del cluster_atoms[cluster_mole], cluster_pos[cluster_mole]
                if cluster_mole in self._cation_moles:
                    recorded_cations.append(cluster_mole)
                    ligand_molecules[cluster_mole] = ligand_molecules[anion_mole]
            
            for mol in ligand_molecules[anion_mole]:
                if mol in self._cation_moles:
                    ligand_molecules[mol] = ligand_molecules[anion_mole]
        
        for cation_mole in self._cation_moles:
            if cation_mole in recorded_cations:
                continue
            cluster_moles = self.walking_list(cation_mole, cluster_graph)
            cen_pos = molecule_position[cation_mole]
            for cluster_mole in cluster_moles:
                ligand_molecules[cation_mole].append(cluster_mole)
                cluster_atoms[cation_mole].extend(molecule_dictionary[cluster_mole])
                if cluster_mole in cluster_graph[cation_mole].keys():
                    for _ in range(cluster_graph[cation_mole][cluster_mole]):
                        ligand_molecules[cation_mole].append(cluster_mole)
                mea_pos = molecule_position[cluster_mole]
                dist_mat = self.distance_matrix_wrap(cen_pos, mea_pos, self._AtomStep.lattice, self._period_images)
                for _, jj, kk in zip(*np.where(dist_mat < self._lattice_bond)):
                    mea_pos[jj] += self._period_images[kk]
                    mea_pos = self.image_shift(mea_pos, self._AtomStep.lattice, benchmark=jj, cartesian=0)
                    break
                cluster_pos[cation_mole].extend(list(mea_pos))
        for mol, positions in cluster_pos.items():
            cluster_pos[mol] = list(self._AtomStep.f2c(array(positions), self.lattice) - self._AtomStep.f2c(np.average(positions, axis=0), self.lattice))
        return ligand_molecules, cluster_atoms, cluster_pos
    def walking_list(self, index, graph):
        new, mole_list = [index], []
        while(new != []):
            atom = new.pop(0)
            if atom not in mole_list:
                mole_list.append(atom)
            if atom in self._solvent_moles:
                continue
            elif atom in self._anion_moles:
                connected_atom = list(atom2 for atom2 in graph[atom].keys() if graph[atom][atom2] and atom2 in self._cation_moles)
            else:
                connected_atom = list(atom2 for atom2 in graph[atom].keys() if graph[atom][atom2])
            new.extend([atom for atom in connected_atom if atom not in mole_list])
        return mole_list

if __name__ == "__main__":
    arg = args({"output": ["ligand", "output file name includes bonded atoms, distances to the center, and bonded molecules"], \
                "output2": ["clusters", "output directory name contains .xyz files, each representing a cluster"], \
                "output3": ["cluster", "output file name includes cluster type, count, percentage, and center atom"], \
                "output4": ["POSCAR_first_shield", "POSCAR file for the first shield"], \
                "mode": ["a", str, "atom mode, molecule mode, or SSIP/CIP/AGG mode"], \
                "bond": [None, str, "specify a file with user-defined bond length threshold"], \
                "molecules": [None, str, "specify a file with molecule information"], \
                "scaling": [2.5, float, "set larger if you have large molecules"]}, "POSCAR")
    if poscar.args.mode.lower() not in ["a", "m", "g"]:
        print("Warning: mode must be \"a\" (atom mode), \"m\" (molecule mode), or \"g\" (SSIP/CIP/AGG mode). Reset as \"a\".")
        poscar.args.mode = "a"

    poscar = POSCAR(arg)
    if poscar.args.atomfile is None:
        print("\nWhich atoms you would like to catch?")
    
    atom_list = line(poscar, atom=Atom())
    poscar.read_all(atoms=atom_list)

    if poscar.args.mode.lower() == "a" or poscar.args.mode.lower() == "m":
        print("\nCenter atom:")
        cen_atom = line(poscar, atom=Atom())
    elif poscar.args.mode == "g":
        cen_atom = Atom(poscar)
        while(1):
            tmp = input("\nCation: ")
            cation = Atom(poscar)
            for _ in tmp.split(","):
                err = cation.put(_)
                if err < 0:
                    print(cation.err_list[err])
                    break
            else:
                cen_atom.put(cation)
                break
        while(1):
            tmp = input("Anion: ")
            anion = Atom(poscar)
            for _ in tmp.split(","):
                err = anion.put(_)
                if err < 0:
                    print(anion.err_list[err])
                    break
            else:
                break
        cen_atom.put(anion)

    center_atom_step, cluster_atom_step = AtomStep_Single_Point(poscar), AtomStep_Single_Point(poscar)
    center_pos, cluster_pos = [], []
    center_atom_step.atoms = cen_atom
    cluster_atom = Atom(poscar)
    if poscar.args.mode == "g":
        for atom, pos in zip(atom_list.get(), poscar.fractional_position):
            if atom in cen_atom.get():
                center_pos.append(pos)
            cluster_pos.append(pos)
            cluster_atom.put(atom)
    else:
        for atom, pos in zip(atom_list.get(), poscar.fractional_position):
            if atom not in cen_atom.get():
                cluster_pos.append(pos)
                cluster_atom.put(atom)
            else:
                center_pos.append(pos)
    cluster_atom_step.atoms = cluster_atom
    center_atom_step.lock(); cluster_atom_step.lock()
    center_atom_step.fractional_position, cluster_atom_step.fractional_position = center_pos, cluster_pos
    if poscar.args.mode == "g":
        cluster = HCE(poscar)
    else:
        cluster = Ligand(poscar)
    while(1):
        tmp = input("Lattice bond [2.5]: ")
        if tmp == "":
            cluster.lattice_bond = 2.5
            break
        elif compile(r"^\d+\.?\d*$").match(str(tmp)) is not None:
            cluster.lattice_bond = float(tmp)
            break
        else:
            print("Warning: Input a positive number.")
    
    cluster.center_atom_step, cluster.measure_atom_step = center_atom_step, cluster_atom_step
    
    if poscar.args.molecules is not None and path.isfile(poscar.args.molecules):
        cluster.read_molecules(poscar.args.molecules)
    elif poscar.args.molecules is not None:
        print(f"Warning: {poscar.args.molecules} does not exist.")
    if poscar.args.bond is not None:
        cluster.molecules.read_bond(poscar.args.bond)
    else:
        bond_type = threshold()
        if bond_type is not None:
            cluster.atom_step.molecules.threshold = bond_type
    
    cluster.build_graph(molecule=True if cluster.molecules.molecule_dictionary != {} else False)
    cluster.graph_to_molecule_position()

    if poscar.args.mode == "g":
        cluster.ligand_search(cation, anion)
        anion_mole, cation_mole = [], []
        for mol, atoms in cluster.molecules.molecule_dictionary.items():
            for atom in atoms:
                if atom in cation.get():
                    cation_mole.append(mol)
                    break
                elif atom in anion.get():
                    anion_mole.append(mol)
                    break
    else:
        cluster.ligand_search()
    
    if not path.isdir(poscar.args.output2):
        mkdir(poscar.args.output2)

    mole_kind = [kind.split('_')[0] for kind in cluster.molecules.molecule_kind]
    mole_kind_sorted = sorted(list(set(mole_kind)))
    
    for cen_mol, atoms in cluster.cluster_atoms.items():
        elements = []
        for atom in atoms:
            elements.append(poscar.elements[poscar.atom_step.index_list[atom]])
        if cluster.molecules.molecule_kind != []:
            mole_kind, cluster_type = cluster.cluster_types()
            filename = "".join([f"{kind}{type}" if i == 0 else f"_{kind}{type}" for i, (kind, type) in enumerate(zip(mole_kind, cluster_type[cen_mol]))])
        else:
            filename = cen_mol
        
        if poscar.args.mode == "g":
            for atom in cluster.cluster_atoms[cen_mol]:
                if atom in cation.get():
                    break
            else:
                continue
        
        with open(f"{poscar.args.output2}/{filename}.xyz", "a") as write_file:
            write_file.write(f"{len(elements)}\n")
            write_file.write("center")
            for atom in cluster.molecules.molecule_dictionary[cen_mol]:
                write_file.write(f",{atom+poscar.args.atom}")
            write_file.write(",,measure")

            for atom in cluster.cluster_atoms[cen_mol][len(cluster.molecules.molecule_dictionary[cen_mol]):]:
                write_file.write(f",{atom+poscar.args.atom}")
            write_file.write("\n")
            for idx, (element, pos) in enumerate(zip(elements, cluster.cluster_position[cen_mol])):
                write_file.write(f"{element}\t{pos[0]:.8f}\t{pos[1]:.8f}\t{pos[2]:.8f}\n")
    
    with open(f"{poscar.args.output}.csv", "w") as write_file:
        write_file.write(f"center atom,element,molecule,bonded atom,element,molecule,distance")
        if poscar.args.mode == "g":
            write_file.write(",cluster type\n")
        else:
            write_file.write("\n")
        
        recorded_cen_moles = []
        if poscar.args.mode == "g":
            hce_cluster_type = {}
        cluster_atoms, cluster_position = cluster.cluster_atoms, cluster.cluster_position
        neighbors = cluster.build_neighbors(cluster.cluster_graph)
        
        for cen_mol in sorted(cluster.cluster_graph.keys()):
            if cen_mol in recorded_cen_moles:
                continue
            recorded_cen_moles.append(cen_mol)
            cen_atoms_idx = cluster.molecules.molecule_dictionary[cen_mol]
            if cluster.args.mode != "g" and (len(cen_atoms_idx) == len(cluster_atoms[cen_mol]) or cluster_atoms[cen_mol] == []):
                continue
            elif cluster.args.mode == "g" and cen_mol not in cation_mole:
                continue
            elif cluster.args.mode == "g" and (len(cen_atoms_idx) == len(cluster_atoms[cen_mol]) or cluster_atoms[cen_mol] == []):
                cluster_atoms[cen_mol].extend(cluster.molecules.molecule_dictionary[cen_mol])
                cluster_position[cen_mol].extend(cluster.molecules.molecule_position[cen_mol])
                recorded_neighbor_moles = []
                for neighbor in neighbors[cen_mol]:
                    if neighbor in recorded_neighbor_moles:
                        continue
                    recorded_neighbor_moles.append(neighbor)
                    cluster_atoms[cen_mol].extend(cluster.molecules.molecule_dictionary[neighbor])
                    mea_pos = cluster.molecules.molecule_position[neighbor]
                    dist_mat = cluster.distance_matrix_wrap(cluster.molecules.molecule_position[cen_mol], mea_pos, cluster.lattice, cluster.period_images)
                    for _, jj, kk in zip(*np.where(dist_mat < cluster.lattice_bond)):
                        mea_pos[jj] += cluster.period_images[kk]
                        mea_pos = cluster.image_shift(mea_pos, cluster.lattice, benchmark=jj, cartesian=0)
                        break
                    cluster_position[cen_mol].extend(mea_pos)
                cluster_position[cen_mol] = cluster.atom_step.f2c(array(cluster_position[cen_mol]), cluster.lattice)
            cen_pos, mea_pos = cluster.atom_step.c2f(array(cluster_position[cen_mol][:len(cen_atoms_idx)]), cluster.lattice), cluster.atom_step.c2f(array(cluster_position[cen_mol][len(cen_atoms_idx):]), cluster.lattice)
            dist_vec = cluster.distance_matrix_func(cen_pos, mea_pos, cluster.lattice)
            for i, j in zip(*where(dist_vec < cluster.lattice_bond)):
                cen_idx, mea_idx = cluster_atoms[cen_mol][i], cluster_atoms[cen_mol][len(cen_atoms_idx)+j]
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

                write_file.write(f"{cen_atoms_idx[i]+poscar.args.atom},{poscar.elements[poscar.atom_step.index_list[cen_atoms_idx[i]]]},")
                if cluster.molecules.molecule_kind != []:
                    write_file.write(f"{cluster.molecules.molecule_kind[cen_mol]},")
                else:
                    write_file.write(f"{cen_mol+1},")
                write_file.write(f"{cluster.cluster_atoms[cen_mol][len(cen_atoms_idx)+j]+poscar.args.atom},{poscar.elements[poscar.atom_step.index_list[cluster.cluster_atoms[cen_mol][len(cen_atoms_idx)+j]]]},")
                atom_idx = poscar.atom_step.index_list[cluster.cluster_atoms[cen_mol][len(cen_atoms_idx)+j]]
                for ligand in cluster.ligand_molecules[cen_mol]:
                    if atom_idx in cluster.molecules.molecule_dictionary[ligand]:
                        if cluster.molecules.molecule_kind != []:
                            write_file.write(f"{cluster.molecules.molecule_kind[ligand]}")
                            break
                        else:
                            write_file.write(f"{ligand+1}")
                            break
                write_file.write(f",{dist_vec[i, j]:.4f}")
                if poscar.args.mode == "g":
                    anion_counter, cation_counter = 0, 0
                    recorded_an_cat_mole = []
                    for mol in cluster.ligand_molecules[cen_mol]:
                        if mol in anion_mole and mol not in recorded_an_cat_mole:
                            anion_counter += 1
                        elif mol in cation_mole and mol not in recorded_an_cat_mole:
                            cation_counter += 1
                        recorded_an_cat_mole.append(mol)
                    if anion_counter == 0:
                        hce_cluster_type[cen_mol] = "SSIP"
                        write_file.write(",SSIP")
                    elif cation_counter == 1:
                        hce_cluster_type[cen_mol] = "CIP"
                        write_file.write(",CIP")
                    elif cation_counter > 1:
                        hce_cluster_type[cen_mol] = "AGG"
                        write_file.write(",AGG")
                write_file.write("\n")

    if cluster.molecules.molecule_kind != []:
        with open(f"{poscar.args.output3}.csv", "w") as write_file:
            cluster_sum = 0
            cluster_type_count = cluster.defaultdict(int)
            for mol, ligand_count in cluster_type.items():
                if ligand_count is None:
                    continue
                cluster_sum += 1
                ligand_count_str = str(ligand_count).replace('[', '').replace(']', '').replace(',', '')
                cluster_type_count[ligand_count_str] += 1
            
            kinds = "".join([str(_) if i == 0 else f",{_}" for i, _ in enumerate(mole_kind)])
            write_file.write(f"{kinds},sum,percentage\n")
            for key, value in cluster_type_count.items():
                write_file.write(f"{key.replace(' ', ',')},{value},{value / cluster_sum:.4f}\n")
            write_file.write('\n')

            if poscar.args.mode == "g":
                from collections import Counter
                hce_cluster_counter = Counter(hce_cluster_type.values())
                atoms_sum = len(hce_cluster_type.values())
                write_file.write(f"SSIP,CIP,AGG\n{hce_cluster_counter['SSIP']},{hce_cluster_counter['CIP']},{hce_cluster_counter['AGG']}\n{hce_cluster_counter['SSIP']/atoms_sum:.4f},{hce_cluster_counter['CIP']/atoms_sum:.4f},{hce_cluster_counter['AGG']/atoms_sum:.4f}\n\n")

            write_file.write(f"center atom,{kinds}\n")
            for key in sorted(cluster_type.keys()):
                if poscar.args.mode == "g" and key not in cation_mole:
                    continue
                mole_dict_str = str(cluster.molecules.molecule_dictionary[key]).replace('[', '').replace(']', '').replace(',', '')
                write_file.write(f"{mole_dict_str}")
                for i in cluster_type[key]:
                    write_file.write(f",{i}")
                write_file.write("\n")

    if poscar.args.mode == "g":
        poscar_first_solvation = POSCAR(arg)
        poscar_first_solvation.args.output = poscar.args.output4
        recorded_moles, positions, elements = [], [], []
        for moles in cluster.ligand_molecules.values():
            for mole in moles:
                if mole not in recorded_moles:
                    recorded_moles.append(mole)
                    for atom in cluster.molecules.molecule_dictionary[mole]:
                        elements.append(poscar.elements[poscar.atom_step.index_list[atom]])
                    positions.extend(cluster.molecules.molecule_position[mole])
        poscar_first_solvation.lattice = poscar.lattice
        poscar_first_solvation.rearrange(elements, positions)
        poscar_first_solvation.write_all()
    print("Done!")
