import numpy as np
from kit.fundamental import Atom, Graph
from kit.interface import args, step, threshold
from kit.function import Solvent
from kit.vasp import XDATCAR

class Reaction(Solvent):
    def __init__(self, input_obj=None):
        Solvent.__init__(self, input_obj)
    @property
    def bond_info(self):
        return self._BondInfo
    def bond_breaking(self, t):
        mole_dict = self._fragment_mole_dict
        for mol, position in self._fragment_mole_pos.items():
            dist_mat = self.distance_matrix_func(position[t], position[t], self.lattice)
            for idx_i, cen_atom in enumerate(mole_dict[mol]):
                cen_atom_idx = self._AtomStep.atoms.index_list[cen_atom]
                for mea_atom_idx in np.where(self._adj_mat[cen_atom_idx] > 0)[0]:
                    mea_atom = self._AtomStep.atom_list[mea_atom_idx]
                    idx_j = mole_dict[mol].index(self._AtomStep.atom_list[mea_atom_idx])
                    if idx_i < idx_j:
                        continue
                    if dist_mat[idx_i, idx_j] > self._adj_mat[cen_atom_idx, mea_atom_idx]:
                        self._adj_mat = self.delete_bond(cen_atom, mea_atom, self._adj_mat)
                        atom_min, atom_max = min(cen_atom, mea_atom), max(cen_atom, mea_atom)
                        fragments = self.connect(cen_atom, mea_atom, self._adj_mat)
                        position_min = self._AtomStep.f2c(position[t][idx_i], self.lattice) if cen_atom == atom_min else self._AtomStep.f2c(position[t][idx_j], self.lattice)
                        position_max = self._AtomStep.f2c(position[t][idx_i], self.lattice) if cen_atom == atom_max else self._AtomStep.f2c(position[t][idx_j], self.lattice)
                        bond_type = 1 if len(fragments) == 1 else 2
                        distance = np.linalg.norm(position_min - position_max)
                        info = [t, atom_min, atom_max, position_min, position_max, distance]
                        if not any([bond_info[0] == info[0] and bond_info[1] == info[1] and bond_info[2] == info[2] for bond_info in self._BondInfo[bond_type]]):
                            self._BondInfo[bond_type].append(info)
        self.update_graph()
    def bond_creating(self, t):
        mole_dict, mole_pos_dict = self._fragment_mole_dict, self._fragment_mole_pos
        recorded = set()
        avg_pos_arr = []
        avg_pos = {}
        for mol, position in mole_pos_dict.items():
            avg_pos[mol] = np.average(position[t], axis=0)
            avg_pos_arr.append(avg_pos[mol])
        avg_pos_arr = np.array(avg_pos_arr)
        dist_mat = self.distance_matrix_wrap(avg_pos_arr, avg_pos_arr, self.lattice, self._period_images)
        avg_pos_keys = list(avg_pos.keys())
        elements = self._AtomStep.atoms.elements
        for mole_idx_i, mole_idx_j, period_image in zip(*np.where((dist_mat < max(self._AtomStep.molecules.threshold.values())*self._search_scaling) & (dist_mat > 0.01))):
            #if mole_idx_i < mole_idx_j:
            #    continue
            mole_i, mole_j = avg_pos_keys[mole_idx_i], avg_pos_keys[mole_idx_j]
            frac_pos_i, frac_pos_j = mole_pos_dict[mole_i][t], mole_pos_dict[mole_j][t] + self._period_images[period_image]
            mole_dist_mat = self.distance_matrix_func(frac_pos_i, frac_pos_j, self.lattice)
            for idx_i, idx_j in zip(*np.where((mole_dist_mat < max(self._AtomStep.molecules.threshold.values())) & (mole_dist_mat > 0.01))):
                if idx_i < idx_j:
                    continue
                cen_atom, mea_atom = mole_dict[mole_i][idx_i], mole_dict[mole_j][idx_j]
                cen_atom_idx, mea_atom_idx = self._AtomStep.atoms.index_list[cen_atom], self._AtomStep.atoms.index_list[mea_atom]
                if 0.01 < mole_dist_mat[idx_i, idx_j] <= self._AtomStep.molecules.threshold[f"{elements[cen_atom_idx]}-{elements[mea_atom_idx]}"]:
                    self._adj_mat = self.create_bond(cen_atom, mea_atom, self._adj_mat)
                    atom_min, atom_max = min(cen_atom, mea_atom), max(cen_atom, mea_atom)
                    position_min = self._AtomStep.f2c(frac_pos_i[idx_i], self.lattice) if cen_atom == atom_min else self._AtomStep.f2c(frac_pos_j[idx_j], self.lattice)
                    position_max = self._AtomStep.f2c(frac_pos_i[idx_i], self.lattice) if cen_atom == atom_max else self._AtomStep.f2c(frac_pos_j[idx_j], self.lattice)
                    distance = np.linalg.norm(position_min - position_max)
                    info = [t, atom_min, atom_max, position_min, position_max, distance]
                    if not any([bond_info[0] == info[0] and bond_info[1] == info[1] and bond_info[2] == info[2] for bond_info in self._BondInfo[3]]):
                        self._BondInfo[3].append(info)
                        if not np.all(self._period_images[period_image] == np.array([0, 0, 0])):
                            recorded.add((mole_j, period_image, idx_j))
        for info in list(recorded):
            mole_pos_dict[info[0]][t][info[2]] += self._period_images[info[1]]
            mole_pos_dict[info[0]][t] = self.image_shift(mole_pos_dict[info[0]][t], self.lattice, benchmark=info[2], cartesian=0)
            mole_pos_dict[info[0]][t:] = self.image_shift(mole_pos_dict[info[0]][t:], self.lattice, cartesian=0)
        for mol, position in zip(mole_dict.values(), mole_pos_dict.values()):
            if len(position[t]) < 2:
                continue
            dist_mat = self.distance_matrix_func(position[t], position[t], self.lattice)
            for idx_i, idx_j in zip(*np.where((dist_mat < max(self._AtomStep.molecules.threshold.values())*self._search_scaling) & (dist_mat > 0.01))):
                if idx_i < idx_j:
                    continue
                cen_atom, mea_atom = mol[idx_i], mol[idx_j]
                cen_atom_idx, mea_atom_idx = self._AtomStep.atoms.index_list[cen_atom], self._AtomStep.atoms.index_list[mea_atom]
                if 0.01 < dist_mat[idx_i, idx_j] <= self._AtomStep.molecules.threshold[f"{elements[cen_atom_idx]}-{elements[mea_atom_idx]}"] and self._adj_mat[cen_atom_idx, mea_atom_idx] == 0:
                    self._adj_mat[cen_atom_idx, mea_atom_idx] = self._adj_mat[mea_atom_idx, cen_atom_idx] = self._AtomStep.molecules.threshold[f"{elements[cen_atom_idx]}-{elements[mea_atom_idx]}"]
                    atom_min, atom_max = min(cen_atom, mea_atom), max(cen_atom, mea_atom)
                    position_min = self._AtomStep.f2c(position[t][idx_i], self.lattice) if cen_atom == atom_min else self._AtomStep.f2c(position[t][idx_j], self.lattice)
                    position_max = self._AtomStep.f2c(position[t][idx_i], self.lattice) if cen_atom == atom_max else self._AtomStep.f2c(position[t][idx_j], self.lattice)
                    distance = np.linalg.norm(position_min - position_max)
                    info = [t, atom_min, atom_max, position_min, position_max, distance]
                    if not any([bond_info[0] == info[0] and bond_info[1] == info[1] and bond_info[2] == info[2] for bond_info in self._BondInfo[3]]):
                        self._BondInfo[3].append(info)
        self.update_graph()
    def run(self):
        if len(self._AtomStep.molecules.molecule_position[0]) > 1000:
            print("Completed step:")
        self._BondInfo = {1: [], 2: [], 3: []}
        for t in range(len(self._AtomStep.molecules.molecule_position[0])):
            if t > 0 and t % 1000 == 0:
                print(f"{t+min(self._AtomStep.steps.get())} / {max(self._AtomStep.steps.get())}")
            self.bond_breaking(t)
            self.bond_creating(t)
        self._BondInfo[1].sort(key=lambda x: x[0]); self._BondInfo[2].sort(key=lambda x: x[0]); self._BondInfo[3].sort(key=lambda x: x[0])

if __name__ == "__main__":
    arg = args({"output": "reaction", "molecules": [None, str, "specify a file with molecule information"], \
                "bond": [None, str, "specify a file with bond threshold"], \
                "scaling": [2.5, float, "set larger if you have large molecules"]}, "XDATCAR")
    
    xdatcar = XDATCAR(arg)
    
    step_list = step(xdatcar)
    if arg.args.bond is not None:
        xdatcar.molecules.read_bond(arg.args.bond)
    else:
        bond_type = threshold()
        if bond_type is not None:
            xdatcar.atom_step.molecules.threshold = bond_type
    
    xdatcar.read_molecules(arg.args.molecules)
    atom_list = Atom(xdatcar)
    for mole_list in xdatcar.molecules.molecule_dictionary.values():
        for atom in mole_list:
            atom_list.put(atom)
    atom_list.sort()
    
    print("Reading XDATCAR...")
    xdatcar.read_all(steps=step_list, atoms=atom_list)
    xdatcar.atom_step.build_molecule_position()
    
    ReactionInit = Graph(xdatcar)
    reaction = Reaction(xdatcar)
    reaction.period_images = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0],
                                       [0, 1, 0], [0, -1, 0], [1, 1, 0],
                                       [1, -1, 0], [-1, 1, 0], [-1, -1, 0]])
    ReactionInit.elements = reaction.elements = [xdatcar.elements[element] for element in atom_list.get()]
    
    ReactionInit.build_graph(wrap=False)
    reaction.adjacent_matrix = ReactionInit.adjacent_matrix

    print("Analysing reaction...")
    reaction.run()

    elements = reaction.atom_step.atoms.elements
    idx_list = reaction.atom_step.atoms.index_list
    with open("bond_breaking.dat", "w") as write_file:
        if reaction.bond_info[1] + reaction.bond_info[2] == []:
            write_file.write("No bond is broken\n")
        else:
            for i, mole_list in reaction.molecules.molecule_dictionary.items():
                counter = 0
                write_file.write(f"Molecule: {reaction.molecules.molecule_kind[i]}\n")
                removes = []
                for j, info in enumerate(reaction.bond_info[1]):
                    if info[1] in mole_list or info[2] in mole_list:
                        if any([bond_info[1] == info[1] and bond_info[2] == info[2] for bond_info in reaction.bond_info[1][j+1:]]):
                            removes.append(info)
                            continue
                        write_file.write(f"\tStep: {info[0]}, New fragments: No\n")
                        write_file.write(f"\t\tAtom {info[1]}, Element: {elements[idx_list[info[1]]]}, Position: ({info[3][0]:.2f}, {info[3][1]:.2f}, {info[3][2]:.2f})\n")
                        write_file.write(f"\t\tAtom {info[2]}, Element: {elements[idx_list[info[2]]]}, Position: ({info[4][0]:.2f}, {info[4][1]:.2f}, {info[4][2]:.2f})\n")
                        write_file.write(f"\t\t{elements[idx_list[info[1]]]}-{elements[idx_list[info[2]]]} threshold: {reaction.molecules.threshold[f'{elements[idx_list[info[1]]]}-{elements[idx_list[info[2]]]}']:.2f} A. Distance: {info[5]:.2f} A\n")
                        counter += 1
                        removes.append(info)
                if removes != []:
                    for remove in removes:
                        reaction.bond_info[1].remove(remove)
                removes = []
                for j, info in enumerate(reaction.bond_info[2]):
                    if info[1] in mole_list or info[2] in mole_list:
                        for build_info in reaction.bond_info[3]:
                            if info[1] == build_info[1] and info[2] == build_info[2] and build_info[0] > info[0]:
                                reaction.bond_info[3].remove(build_info)
                                break
                        else:
                            if any([bond_info[1] == info[1] and bond_info[2] == info[2] for bond_info in reaction.bond_info[2][j+1:]]):
                                removes.append(info)
                                continue
                            write_file.write(f"\tStep: {info[0]}, New fragments: Yes\n")
                            write_file.write(f"\t\tAtom {info[1]}, Element: {elements[idx_list[info[1]]]}, Position: ({info[3][0]:.2f}, {info[3][1]:.2f}, {info[3][2]:.2f})\n")
                            write_file.write(f"\t\tAtom {info[2]}, Element: {elements[idx_list[info[2]]]}, Position: ({info[4][0]:.2f}, {info[4][1]:.2f}, {info[4][2]:.2f})\n")
                            write_file.write(f"\t\t{elements[idx_list[info[1]]]}-{elements[idx_list[info[2]]]} threshold: {reaction.molecules.threshold[f'{elements[idx_list[info[1]]]}-{elements[idx_list[info[2]]]}']:.2f} A. Distance: {info[5]:.2f} A\n")
                            counter += 1
                        removes.append(info)
                if removes != []:
                    for remove in removes:
                        reaction.bond_info[2].remove(remove)
                write_file.write(f"Break {counter} bond")
                if counter > 1:
                    write_file.write("s")
                write_file.write("\n\n")
    with open("bond_building.dat", "w") as write_file:
        if reaction.bond_info[3] == []:
            write_file.write("No bond is built\n")
        else:
            for i, mole_list in reaction.molecules.molecule_dictionary.items():
                counter = 0
                write_file.write(f"Molecule: {reaction.molecules.molecule_kind[i]}\n")
                removes = []
                for j, info in enumerate(reaction.bond_info[3]):
                    if info[1] in mole_list or info[2] in mole_list:
                        if any([bond_info[1] == info[1] and bond_info[2] == info[2] for bond_info in reaction.bond_info[3][j+1:]]):
                            removes.append(info)
                            continue
                        for j, inner_mole_list in reaction.molecules.molecule_dictionary.items():
                            if info[1] in inner_mole_list:
                                MoleMin = reaction.molecules.molecule_kind[j]
                            if info[2] in inner_mole_list:
                                MoleMax = reaction.molecules.molecule_kind[j]
                        write_file.write(f"\tStep: {info[0]}\n")
                        write_file.write(f"\t\tAtom {info[1]}, Element: {elements[idx_list[info[1]]]}, Position: ({info[3][0]:.2f}, {info[3][1]:.2f}, {info[3][2]:.2f}), Molecule: {MoleMin}\n")
                        write_file.write(f"\t\tAtom {info[2]}, Element: {elements[idx_list[info[2]]]}, Position: ({info[4][0]:.2f}, {info[4][1]:.2f}, {info[4][2]:.2f}), Molecule: {MoleMax}\n")
                        write_file.write(f"\t\t{elements[idx_list[info[1]]]}-{elements[idx_list[info[2]]]} threshold: {reaction.molecules.threshold[f'{elements[idx_list[info[1]]]}-{elements[idx_list[info[2]]]}']:.2f} A. Distance: {info[5]:.2f} A\n")
                        counter += 1
                        removes.append(info)
                if removes != []:
                    for remove in removes:
                        reaction.bond_info[3].remove(remove)
                write_file.write(f"Build {counter} bond")
                if counter > 1:
                    write_file.write("s")
                write_file.write("\n\n")
    print("Done!")
