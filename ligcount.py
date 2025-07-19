from os import getcwd, cpu_count, mkdir, path
from collections import defaultdict
from re import compile
from concurrent.futures import ProcessPoolExecutor
from collections import defaultdict
from numpy import array, array_split
from kit.fundamental import AtomStep_Trj
from kit.function import Ligand
from kit.interface import args, blocks, step_lines
from kit.vasp import XDATCAR
from ligand import HCE

def function(num, total_AtomStep, center_AtomStep, measure_AtomStep, cutoff, electrolyte):
    if center_AtomStep.args.mode.lower() == "g":
        molecule = HCE()
    else:
        molecule = Ligand()
    molecule.args = center_AtomStep.args
    molecule.atom_step, molecule.center_atom_step, molecule.measure_atom_step = total_AtomStep, center_AtomStep, measure_AtomStep
    molecule.lattice_bond = cutoff
    molecule.search_scaling = total_AtomStep.args.scaling
    print(f"Analysing {num}...")
    molecule.read_molecules(total_AtomStep.args.molecules)
    molecule.atom_step.build_molecule_position()
    molecule.ligand_search(cation=electrolyte[0], anion=electrolyte[1])
    return center_AtomStep.args.mode.lower(), molecule

if __name__ == "__main__":
    arg = args({"output": "count", "output2": ["clusters_trj", "output directory name contains .xyz files, each representing a cluster"], \
                "molecules": [None, str, "specify a file with molecule information"], "scaling": [2.5, float, "set larger if you have large molecules"], \
                "mode": ["a", str, "atom mode, molecule mode, or SSIP/CIP/AGG mode"], \
                "cpu": [int(cpu_count()/2.5), int, "specify the maximum number of applied CPU cores"]}, "XDATCAR")
    
    xdatcar = XDATCAR(arg)
    
    atom_list, atom_list_flatten, annotate_info = blocks(xdatcar, 2, flatten_flag=True, annotate=1)
    flag = False
    for mode in annotate_info:
        if mode == []:
            mode.append(arg.args.mode.lower())
        else:
            flag = True
    if arg.args.mode == "g":
        flag = True
    step_list, step_list_flatten = step_lines(xdatcar, items=len(atom_list), flatten_flag=True)
    while(arg.args.molecules is None or not path.isfile(arg.args.molecules)):
        print("\nWarning: 'molecules' argument is not a file.")
        tmp = input("Specify a file with molecule information: ")
        arg.args.molecules = tmp
    
    if flag:
        cation_mole, anion_mole = [], []
        from kit.fundamental import Atom
        cen_atom = Atom(xdatcar)
        while(1):
            tmp = input("Cation: ")
            cation = Atom(xdatcar)
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
            anion = Atom(xdatcar)
            for _ in tmp.split(","):
                err = anion.put(_)
                if err < 0:
                    print(anion.err_list[err])
                    break
            else:
                cen_atom.put(anion)
                break
        
        cluster_atom = Atom(xdatcar)
        for atom in atom_list_flatten.get():
            cluster_atom.put(atom)
        
        for mol, atoms in xdatcar.molecules.molecule_dictionary.items():
            for atom in atoms:
                if atom in cation.get():
                    cation_mole.append(mol)
                    break
                elif atom in anion.get():
                    anion_mole.append(mol)
                    break
    while(1):
        tmp = input("\nCutoff radius [2.5]: ")
        if tmp == '':
            cutoff = 2.5
            break
        elif compile(r"\d+\.?\d*").match(tmp) is not None:
            cutoff = float(tmp)
            break
        else:
            print("Warning: Input error.")
    while(1):
        plot_flag = False
        tmp = input("Plot (y/n) [n]: ")
        if tmp.lower() == '':
            plot_flag = False
            break
        elif tmp.lower() == 'y':
            plot_flag = True
            import matplotlib.pyplot as plt
            from seaborn import color_palette
            break
        elif tmp.lower() == 'n':
            plot_flag = False
            break
        else:
            print("Warning: Input error.")
    
    print("\nReading XDATCAR...")
    xdatcar.read_all(steps=step_list_flatten, atoms=atom_list_flatten)
    atom_list_flatten_ = atom_list_flatten.get()
    step_list_flatten = (step_list_flatten.get() if step_list_flatten != "all" else "all")
    process_per_line = (int(arg.args.cpu/len(step_list)) if (step_list_flatten == "all" or len(step_list_flatten)*2.5/cpu_count()>15.0) else 1)
    nums = []
    center_AtomSteps, measure_AtomSteps = [], []
    for idx, atoms, steps, mode in zip(range(len(step_list)), atom_list, step_list, annotate_info):
        if mode[0] == "g":
            atoms_1, atoms_2 = cen_atom.get(), cluster_atom.get()
        else:
            atoms_1, atoms_2 = atoms[0].get(), atoms[1].get()
        steps_ = steps.get()
        tmp = []
        for i in steps_:
            tmp.append((slice(0, len(list(range(i.start, i.stop, i.step if i.step is not None else 1))))) if isinstance(i, slice) else i)
        steps_ = tuple(tmp)
        cen_pos, mea_pos = [], []
        for i, atomIdx in enumerate(atom_list_flatten_):
            if atomIdx in atoms_1:
                cen_pos.append(array(xdatcar.fractional_position[:, i]))
            if atomIdx in atoms_2:
                mea_pos.append(array(xdatcar.fractional_position[:, i]))
        steps_split = array_split(steps.get(slice_flatten=True), process_per_line)
        cen_pos = array_split(array(cen_pos).transpose((1, 0, 2)), process_per_line)
        mea_pos = array_split(array(mea_pos).transpose((1, 0, 2)), process_per_line)
        positions = array_split(array(xdatcar.fractional_position), process_per_line)
        print(153, cen_pos[0].shape, len(atoms_1))
        for i in range(process_per_line):
            print(f"Assigning steps {idx+1}_{i+1}...")
            nums.append(f"{idx+1}_{i+1}")
            center_AtomStep, measure_AtomStep = AtomStep_Trj(xdatcar), AtomStep_Trj(xdatcar)
            if mode[0] == "g":
                center_AtomStep.atoms, measure_AtomStep.atoms = cen_atom, cluster_atom
                print(160, len(cen_atom.get()))
            else:
                center_AtomStep.atoms, measure_AtomStep.atoms = atoms[0], atoms[1]
            center_AtomStep.steps.put(steps_split[i]); measure_AtomStep.steps.put(steps_split[i])
            center_AtomStep.fractional_position, measure_AtomStep.fractional_position = cen_pos[i], mea_pos[i]
            center_AtomStep.args = xdatcar.args
            center_AtomStep.args.mode = annotate_info[idx][0]
            center_AtomSteps.append(center_AtomStep); measure_AtomSteps.append(measure_AtomStep)
    
    with ProcessPoolExecutor(arg.args.cpu) as executor:
        futures = [executor.submit(function, num, xdatcar.atom_step, center_AtomStep, measure_AtomStep, cutoff, [cation, anion] if center_AtomStep.args.mode.lower() == "g" else [None, None]) for num, center_AtomStep, measure_AtomStep in zip(nums, center_AtomSteps, measure_AtomSteps)]
        results = [future.result() for future in futures]

    if not path.isdir(arg.args.output2):
        mkdir(arg.args.output2)
    
    cluster_atoms, cluster_pos, bond_num = [defaultdict(list) for _ in range(len(step_list))], [defaultdict(list) for _ in range(len(step_list))], [defaultdict(int) for _ in range(len(step_list))]
    cluster_type_count, bond_num_count = defaultdict(int), defaultdict(int)
    cluster_sum, line_ = 0, 0
    if results[0][0] == "g":
        hce_cluster_type = {}
    for num, result in zip(nums, results):
        if result[0] == "g":
            hce_cluster_type = {}
        mole_kind_sorted, ligand_types = result[1].cluster_types()
        for mol, ligand_counts in ligand_types.items():
            for ligand_count in ligand_counts:
                if result[0] != "g" or (result[0] == "g" and mol in cation_mole):
                    cluster_sum += 1
                    ligand_count_str = str(ligand_count).replace('[', '').replace(']', '').replace(',', '')
                    cluster_type_count[ligand_count_str] += 1
                    bond_num_count[sum(ligand_count)] += 1
        
        if result[0] == "g":
            counter = 0
            for cen_mol, cluster_moles in result[1].ligand_molecules.items():
                for idx, cluster_mole in enumerate(cluster_moles):
                    if cen_mol in anion_mole:
                        continue
                    recorded = []
                    anion_counter, cation_counter = 0, 0
                    for mol in set(cluster_mole):
                        if mol in anion_mole and mol not in recorded:
                            anion_counter += 1
                        elif mol in cation_mole and mol not in recorded:
                            cation_counter += 1
                        recorded.append(mol)
                    if anion_counter == 0:
                        hce_cluster_type[counter] = "SSIP"
                    elif cation_counter == 1:
                        hce_cluster_type[counter] = "CIP"
                    elif cation_counter > 1:
                        hce_cluster_type[counter] = "AGG"
                    counter += 1
        
        line = int(num.split("_")[0])
        if line_ != line:
            if result[0] == "g":
                from collections import Counter
                hce_cluster_counter = Counter(hce_cluster_type.values())
                atoms_sum = len(hce_cluster_type.values())
            with open(arg.args.output+".csv", "a") as write_file:
                write_file.write(f"line {line}\n")
                mole_kind_sorted_str = str(mole_kind_sorted).replace('[', '').replace(']', '').replace(' ', '').replace("'", "")
                write_file.write(mole_kind_sorted_str+",count,ratio\n")
                for key, val in cluster_type_count.items():
                    write_file.write(f"{key.replace(' ', ',')},{val},{val/cluster_sum:.4f}\n")
                write_file.write("\n")
                write_file.write("ligand number,count,ratio\n")
                for key in sorted(bond_num_count.keys()):
                    write_file.write(f"{key},{bond_num_count[key]},{bond_num_count[key]/cluster_sum:.4f}\n")
                if plot_flag:
                    arg.same_name(getcwd(), f"{line}.png")
                    palette = color_palette("Set2", 8)
                    if len(bond_num_count) > 8:
                        palette += color_palette("husl", len(bond_num_count)-8)
                    plt.pie(list(bond_num_count.values()), labels=list(bond_num_count.keys()), radius=1.2, colors=palette)
                    plt.savefig(f"{i+1}.png", dpi=600, format="png")
                    plt.clf()
                write_file.write("\n")
                if result[0] == "g":
                    write_file.write(f"SSIP,CIP,AGG\n{hce_cluster_counter['SSIP']},{hce_cluster_counter['CIP']},{hce_cluster_counter['AGG']}\n{hce_cluster_counter['SSIP']/atoms_sum:.4f},{hce_cluster_counter['CIP']/atoms_sum:.4f},{hce_cluster_counter['AGG']/atoms_sum:.4f}\n\n")
            cluster_type_count, bond_num_count = defaultdict(int), defaultdict(int)
            cluster_sum = 0
            line_ = line

        for mol in ligand_types.keys():
            for cluster_type, cluster_atoms, cluster_pos in zip(ligand_types[mol], result[1].cluster_atoms[mol], result[1].cluster_position[mol]):
                if not len(cluster_atoms):
                    continue
                center_mole = xdatcar.molecules.molecule_dictionary[mol]
                if len(center_mole) == len(cluster_atoms):
                    continue
                kind = "".join([f"_{num}{mole_kind}" for num, mole_kind in zip(cluster_type, mole_kind_sorted)])
                with open(f"{arg.args.output2}/{line}{kind}.xyz", "a") as write_file:
                    c, m = "".join([f",{atom+arg.args.atom}" for atom in cluster_atoms[:len(center_mole)]]), "".join([f",{atom+arg.args.atom}" for atom in cluster_atoms[len(center_mole):]])
                    write_file.write(f"{len(cluster_atoms)}\ncenter{c},,measure{m}\n")
                    elements = [xdatcar.elements[atom] for atom in cluster_atoms]
                    for element, position in zip(elements, cluster_pos):
                        write_file.write(f"{element}\t{position[0]:.8f}\t{position[1]:.8f}\t{position[2]:.8f}\n")
    
    print("Done!")
