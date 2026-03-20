from kit.args import args
arg = args({"output": "count", "mode": ["a", str, "atom mode (a), molecule mode (m), or cluster mode (c)"], \
            "cpu": [-1, int, "specify the maximum number of applied CPU cores"], "samename": [1, int, "0: without same file name check, 1: with same file name check"], \
            "bond": [None, str, "specify a file with user-defined bond length threshold"], "dative": [None, float, "cutoff distance for dative bond"], \
            "scaling": [2.5, float, "set larger if you have large molecules"], \
            "plot": [None, int, "0: no plot, 1: plot"]}, "XDATCAR", molecules=True, output_isdir=True)
from os import mkdir, path, cpu_count
from collections import defaultdict, Counter
from re import compile
from concurrent.futures import ProcessPoolExecutor
from collections import defaultdict
from numpy import array, array_split
from kit.fundamental import Atom, AtomStep_Trj
from kit.function import Ligand, HCE
from kit.interface import blocks, step_lines, read_bond
from kit.vasp import XDATCAR

def function(num, total_AtomStep, center_AtomStep, measure_AtomStep, electrolyte, electrolyte_moles):
    mode = center_AtomStep.args.mode
    if mode == "c":
        ligand = HCE()
    else:
        ligand = Ligand()
    
    ligand.args = center_AtomStep.args
    ligand.atom_step, ligand.center_atom_step, ligand.measure_atom_step = total_AtomStep, center_AtomStep, measure_AtomStep
    ligand.dative_bond = total_AtomStep.args.dative
    ligand.search_scaling = total_AtomStep.args.scaling
    print(f"Analysing {num}...")
    ligand.read_molecules(total_AtomStep.args.molecules)
    ligand.atom_step.build_molecule_position()
    ligand.run(cation=electrolyte[0], anion=electrolyte[1])
    cluster_type_count, bond_num_count = defaultdict(int), defaultdict(int)
    hce_cluster_type = {}
    ligand_types = ligand.cluster_types()[1]
    for mol, ligand_counts in ligand_types.items():
        for ligand_count in ligand_counts:
            if mode != "c" or (mode == "c" and mol in electrolyte_moles[0]):
                ligand_count_str = " ".join(str(_) for _ in ligand_count)
                cluster_type_count[ligand_count_str] += 1
                bond_num_count[sum(ligand_count)] += 1
    
    if mode == "c":
        counter = 0
        for cen_mol, cluster_moles in ligand.ligand_molecules.items():
            for cluster_mole in cluster_moles:
                if cen_mol in electrolyte_moles[1]:
                    continue
                recorded = []
                anion_cation_counter = {"anion": 0, "cation": 0}
                for mol in set(cluster_mole):
                    if mol in electrolyte_moles[1] and mol not in recorded:
                        anion_cation_counter["anion"] += 1
                    elif mol in electrolyte_moles[0] and mol not in recorded:
                        anion_cation_counter["cation"] += 1
                    recorded.append(mol)
                if anion_cation_counter["anion"] == 0:
                    hce_cluster_type[counter] = "SSIP"
                elif anion_cation_counter["cation"] == 1:
                    hce_cluster_type[counter] = "CIP"
                elif anion_cation_counter["cation"] > 1:
                    hce_cluster_type[counter] = "AGG"
                counter += 1
        hce_cluster_counter = Counter(hce_cluster_type.values())
    return mode, ligand, cluster_type_count, bond_num_count, hce_cluster_counter if mode == "c" else None

if __name__ == "__main__":
    if arg.args.mode.lower() not in ["a", "m", "c"]:
        print("Warning: mode must be \"a\" (atom mode), \"m\" (molecule mode), or \"g\" (SSIP/CIP/AGG mode). Reset as \"a\".")
        arg.args.mode = "a"
    elif not arg.args.mode.islower():
        arg.args.mode = arg.args.mode.lower()
    
    xdatcar = XDATCAR(arg)
    
    atom_list, atom_flatten, atom_info, annotate_info = blocks(xdatcar, 3, flatten_flag=True, print_output_flag=True, annotate=1)

    for idx, mode in enumerate(annotate_info):
        if mode == []:
            mode.append(xdatcar.args.mode)
        elif mode[0].lower() == "c":
            if not mode[0].islower():
                mode[0] = mode[0].lower()
        elif mode[0].lower() not in ["a", "m", "c"]:
            print("Line ", idx+1, ": Warning: mode must be \"a\" (atom mode), \"m\" (molecule mode), or \"g\" (SSIP/CIP/AGG mode). Reset as \"a\".")
            mode[0] = "a"
    
    step_list, step_flatten, step_info = step_lines(xdatcar, items=len(atom_list), flatten_flag=True, print_output_flag=True)

    if xdatcar.args.cpu <= 0:
        xdatcar.args.cpu = 3 if cpu_count() > 4 else 1
    elif xdatcar.args.cpu > cpu_count():
        xdatcar.args.cpu = cpu_count()

    while(xdatcar.args.molecules is None or not path.isfile(xdatcar.args.molecules)):
        print("\nWarning: 'molecules' argument is not a file.")
        inp = input("Specify a file with molecule information: ")
        xdatcar.args.molecules = inp
    
    cation_anion_moles, cation_moles, anion_moles = [], [], []
    cation_anion_atom_list, cluster_atom_list = [], []
    for idx, mode in enumerate(annotate_info):
        if mode[0].lower() == "c":
            cen_atoms = Atom(xdatcar); cluster_atoms = Atom(xdatcar)
            cen_atoms.put(atom_list[idx][0]); cen_atoms.put(atom_list[idx][1])
            cluster_atoms.put(atom_list[idx][0]); cluster_atoms.put(atom_list[idx][1]); cluster_atoms.put(atom_list[idx][2])
            cluster_atom_list.append([cen_atoms, cluster_atoms])
            cation_anion_atom_list.append([atom_list[idx][0], atom_list[idx][1]])
            cation_mole, anion_mole = [], []
            for mol, atoms in xdatcar.molecules.molecule_dictionary.items():
                for atom in atoms:
                    if atom in atom_list[idx][0].get():
                        cation_mole.append(mol)
                        break
                    elif atom in atom_list[idx][1].get():
                        anion_mole.append(mol)
                        break
            cation_moles.append(cation_mole); anion_moles.append(anion_mole)
            cation_anion_moles.append([cation_mole, anion_mole])
        else:
            cluster_atom_list.append([])
            cation_anion_atom_list.append([None, None])
            cation_moles.append([]); anion_moles.append([])
            cation_anion_moles.append([None, None])
        
    if xdatcar.args.dative is None or compile(r"^\d+\.?\d*$").match(str(xdatcar.args.dative)) is None:
        while(1):
            inp = input("\nDative bond length [2.5]: ")
            if inp == '':
                xdatcar.args.dative = 2.5
                break
            elif compile(r"\d+\.?\d*").match(inp) is not None:
                xdatcar.args.dative = float(inp)
                break
            else:
                print("Warning: Input error.")
    
    read_bond(xdatcar)

    while(1):
        plot_flag = False
        if xdatcar.args.plot is None or xdatcar.args.plot not in [0, 1]:
            inp = input("Plot (y/n) [n]: ")
        elif xdatcar.args.plot == 1:
            print("Plot (y/n) [y]: y")
            inp = "y"
        elif xdatcar.args.plot == 0:
            print("Plot (y/n) [y]: n")
            inp = "n"
        if inp.lower() == '':
            plot_flag = False
            break
        elif inp.lower() == 'y':
            plot_flag = True
            import matplotlib.pyplot as plt
            from seaborn import color_palette
            break
        elif inp.lower() == 'n':
            plot_flag = False
            break
        else:
            print("Warning: Input error.")
    
    print("\nReading XDATCAR...")
    xdatcar.read_all(steps=step_flatten, atoms=atom_flatten)
    atom_flatten_list = atom_flatten.get()
    step_flatten_list = (step_flatten.get() if step_flatten != "all" else "all")

    if xdatcar.args.cpu < len(step_list)+1:
        cores = [1 for _ in range(len(step_list))]
    elif step_flatten_list == "all" or len(step_flatten_list)*2.5/xdatcar.args.cpu > 15.0:
        weights = [1 if annotate_info[i][0].lower() != "c" else 2 for i in range(len(step_list))]
        weights /= sum(array(weights))
        cores = [int(xdatcar.args.cpu*weights[i]) for i in range(len(step_list))]
        if sum(cores) > xdatcar.args.cpu:
            while(sum(cores) != xdatcar.args.cpu):
                cores[cores.index(max(cores))] -= 1
        elif sum(cores) < xdatcar.args.cpu:
            while(sum(cores) != xdatcar.args.cpu):
                cores[cores.index(min(cores))] += 1
    else:
        cores = [1 for _ in range(len(step_list))]
    nums = []
    total_AtomSteps, center_AtomSteps, measure_AtomSteps = [], [], []
    cation_anion_atom_list_split, cation_anion_moles_split = [], []
    
    for idx, atoms, steps, mode in zip(range(len(step_list)), atom_list, step_list, annotate_info):
        if mode[0] == "c":
            atoms_1, atoms_2 = cluster_atom_list[idx][0], cluster_atom_list[idx][1]
            atoms_1_tuple, atoms_2_tuple = atoms_1.get(), atoms_2.get()
        else:
            atoms_1, atoms_2 = atoms[0], atoms[1]
            atoms_1_tuple, atoms_2_tuple = atoms_1.get(), atoms_2.get()
        total_AtomStep, center_AtomStep, measure_AtomStep = AtomStep_Trj(xdatcar), AtomStep_Trj(xdatcar), AtomStep_Trj(xdatcar)
        total_AtomStep.split(xdatcar.atom_step, steps, xdatcar.atoms); center_AtomStep.split(xdatcar.atom_step, steps, atoms_1); measure_AtomStep.split(xdatcar.atom_step, steps, atoms_2)
        cen_pos= array(center_AtomStep.fractional_position); mea_pos = array(measure_AtomStep.fractional_position)
        steps_split = array_split(steps.get(slice_flatten=True), cores[idx])
        cen_pos = array_split(cen_pos, cores[idx])
        mea_pos = array_split(mea_pos, cores[idx])
        positions = array_split(total_AtomStep.fractional_position, cores[idx])
        for i in range(cores[idx]):
            print(f"Assigning steps {idx+1}_{i+1}...")
            nums.append(f"{idx+1}_{i+1}")
            total_AtomStep, center_AtomStep, measure_AtomStep = AtomStep_Trj(xdatcar), AtomStep_Trj(xdatcar), AtomStep_Trj(xdatcar)
            center_AtomStep.atoms, measure_AtomStep.atoms = atoms_1, atoms_2
            total_AtomStep.atoms = xdatcar.atoms
            total_AtomStep.steps.put(steps_split[i]); center_AtomStep.steps.put(steps_split[i]); measure_AtomStep.steps.put(steps_split[i])
            total_AtomStep.fractional_position, center_AtomStep.fractional_position, measure_AtomStep.fractional_position = positions[i], cen_pos[i], mea_pos[i]
            total_AtomStep.args = xdatcar.args
            center_AtomStep.args.mode = mode[0]
            total_AtomSteps.append(total_AtomStep); center_AtomSteps.append(center_AtomStep); measure_AtomSteps.append(measure_AtomStep)
            cation_anion_atom_list_split.append(cation_anion_atom_list[idx])
            cation_anion_moles_split.append(cation_anion_moles[idx])
    with ProcessPoolExecutor(max_workers=xdatcar.args.cpu) as executor:
        futures = [executor.submit(function, num, total_AtomStep, center_AtomStep, measure_AtomStep, cation_anion_atoms, moles) for num, total_AtomStep, center_AtomStep, measure_AtomStep, cation_anion_atoms, moles in zip(nums, total_AtomSteps, center_AtomSteps, measure_AtomSteps, cation_anion_atom_list_split, cation_anion_moles_split)]
        results = [future.result() for future in futures]

    if path.isdir(arg.args.output) and arg.args.samename:
        rmdir(arg.args.ouput)
    mkdir(xdatcar.args.output)
    mkdir(path.join(xdatcar.args.output, "clusters_trj"))
    if plot_flag:
        mkdir(path.join(xdatcar.args.output, "fig"))
    
    num_count = 0
    for num, result in zip(nums, results):
        line_present = int(num.split("_")[0])
        num_count += 1
        if num_count == 1:
            cluster_type_count = result[2]
            bond_num_count = result[3]
            hce_cluster_counter = result[4]
        else:
            for key, val in result[2].items():
                cluster_type_count[key] += val
            for key, val in result[3].items():
                bond_num_count[key] += val
            if result[0] == "c":
                for key, val in result[4].items():
                    hce_cluster_counter[key] += val
        
        mole_kind_sorted, ligand_types = result[1].cluster_types()
        if cores[int(num.split("_")[0])-1] == num_count:
            if result[0] == "c":
                atoms_sum = sum(hce_cluster_counter.values())
            cluster_sum = sum(cluster_type_count.values())
            with open(path.join(xdatcar.args.output, "count.csv"), "a") as write_file:
                write_file.write(f"line {line_present}\n")
                if result[0] == "c":
                    atoms_indices = f"atoms,{atom_info[line_present-1][0].replace(',', ';')} {atom_info[line_present-1][1].replace(',', ';')} {atom_info[line_present-1][2].replace(',', ';')}\n"
                else:
                    atoms_indices = f"atoms,{atom_info[line_present-1][0].replace(',', ';')} {atom_info[line_present-1][1].replace(',', ';')}\n"
                step_indices = "steps,"
                if str(step_info[line_present-1]).lower() != "all":
                    for idx, step in enumerate(step_info[line_present-1]):
                        if idx > 0:
                            step_indices += f"_{step.replace(',', ';')}"
                        else:
                            step_indices += f"{step.replace(',', ';')}"
                else:
                    step_indices += step_info[line_present-1].lower()
                write_file.write("mode,")
                if result[0].lower() == "a":
                    write_file.write("atom\n")
                elif result[0].lower() == "m":
                    write_file.write("molecule\n")
                elif result[0].lower() == "c":
                    write_file.write("cluster\n")
                write_file.write(atoms_indices + step_indices + '\n\n')
                mole_kind_sorted_str = ",".join(mole_kind for mole_kind in mole_kind_sorted)
                write_file.write(mole_kind_sorted_str+",count,ratio\n")
                for key, val in cluster_type_count.items():
                    write_file.write(f"{key.replace(' ', ',')},{val},{val/cluster_sum:.4f}\n")
                write_file.write("\n")
                if result[0] != "c":
                    bond_num_count = {key: bond_num_count[key] for key in sorted(bond_num_count.keys())}
                    write_file.write("ligand number,count,ratio\n")
                    for bond_num, count in bond_num_count.items():
                        write_file.write(f"{bond_num-1},{count},{count/cluster_sum:.4f}\n")
                    if plot_flag:
                        palette = color_palette("Set2", 8)
                        if len(bond_num_count) > 8:
                            palette += color_palette("husl", len(bond_num_count)-8)
                        plt.pie(list(bond_num_count.values()), labels=list(bond_num_count.keys()), radius=1.2, colors=palette)
                        plt.savefig(path.join(xdatcar.args.output, "fig", f"{line_present}.png"), dpi=600, format="png")
                        plt.clf()
                else:
                    write_file.write(f"SSIP,CIP,AGG\n{hce_cluster_counter['SSIP']},{hce_cluster_counter['CIP']},{hce_cluster_counter['AGG']}\n{hce_cluster_counter['SSIP']/atoms_sum:.4f},{hce_cluster_counter['CIP']/atoms_sum:.4f},{hce_cluster_counter['AGG']/atoms_sum:.4f}\n\n")
                    if plot_flag:
                        plt.pie(list(hce_cluster_counter.values()), labels=list(hce_cluster_counter.keys()), radius=1.2, colors=palette)
                        plt.savefig(path.join(xdatcar.args.output, "fig", f"{line_present}.png"), dpi=600, format="png")
                        plt.clf()
            num_count = 0
        
        for mol in ligand_types.keys():
            for cluster_type, cluster_atoms, cluster_pos in zip(ligand_types[mol], result[1].cluster_atoms[mol], result[1].cluster_position[mol]):
                if not len(cluster_atoms):
                    continue
                center_mole_atoms = xdatcar.molecules.molecule_dictionary[mol]
                if len(center_mole_atoms) == len(cluster_atoms):
                    continue
                kind = "".join([f"_{num}{mole_kind}" for num, mole_kind in zip(cluster_type, mole_kind_sorted)])
                with open(path.join(xdatcar.args.output, "clusters_trj", f"{line_present}{kind}.xyz"), "a") as write_file:
                    c = "".join([f",{atom+xdatcar.args.atom}" for atom in cluster_atoms[:len(center_mole_atoms)]])
                    m = "".join([f",{atom+xdatcar.args.atom}" for atom in cluster_atoms[len(center_mole_atoms):]])
                    write_file.write(f"{len(cluster_atoms)}\ncenter{c},,measure{m}\n")
                    elements = [xdatcar.elements[atom] for atom in cluster_atoms]
                    for element, position in zip(elements, cluster_pos):
                        write_file.write(f"{element}\t{position[0]:.8f}\t{position[1]:.8f}\t{position[2]:.8f}\n")
        
    print("Done!")
