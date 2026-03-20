from os import path, getcwd, listdir
from kit.args import Args
arg = Args({"input": getcwd(), "output": "XDATCAR", "mode": ["None", str, "the input file mode. c: CONQUEST, g: Gromacs, l: LAMMPS, q: Quantum ESPRESSO"]})
from re import compile
from collections import defaultdict
from numpy import zeros, array
from kit.fundamental import AtomStep_Trj
from kit.software import Conquest, QE
from kit.vasp import XDATCAR

def rearange(atom_info):
    index_dict, ele_num = {}, {}
    elements = []
    
    for val in atom_info.values():
        if val[0] not in elements:
            elements.append(val[0])
    idx = 0
    for element in elements:
        for key, val in atom_info.items():
            if element == val[0]:
                index_dict[idx] = key
                ele_num[element] = 1 if element not in ele_num.keys() else ele_num[element] + 1
                idx += 1
    return ele_num, index_dict

if __name__ == "__main__":
    arg.same_name(getcwd(), "XDATCAR")
    if arg.args.input is None:
        arg.args.input = ""
    if arg.args.mode.lower() not in ["c", "g", "l", "q"] or (not path.isfile(arg.args.input) and not path.isdir(arg.args.input)):
        for f in listdir(getcwd()):
            if "mdtrj" in f:
                arg.args.input = f
                if arg.args.elements is None and path.isfile("QE.in"):
                    arg.args.elements = "QE.in"
                break
            elif "xsf" in f:
                arg.args.input = f
                if arg.args.elements is None and path.isfile("Conquest_input"):
                    arg.args.elements = "Conquest_input"
                break
            elif "gro" in f or "lammpstrj" in f:
                arg.args.input = f
                break
        else:
            arg.args.input = arg.input_file_check(getcwd(), arg.args.input)
    elif path.isfile(arg.args.input):
        arg.args.input = arg.args.input
    elif path.isdir(arg.args.input):
        for f in listdir(arg.args.input):
            if "mdtrj" in f:
                arg.args.input = path.join(arg.args.input, f)
                if arg.args.elements is None and path.isfile(arg.args.input, "QE.in"):
                    arg.args.elements = path.join(arg.args.input, "QE.in")
                break
            elif "xsf" in f:
                arg.args.input = path.join(arg.args.input, f)
                if arg.args.elements is None and path.isfile(arg.args.input, "Conquest_input"):
                    arg.args.elements = path.join(arg.args.input, "Conquest_input")
                break
            elif "gro" in f or "lammpstrj" in f:
                arg.args.input = path.join(arg.args.input, f)
                break
        else:
            arg.args.input = arg.input_file_check(arg.args.input, "")
    else:
        arg.args.input = arg.input_file_check(getcwd(), arg.args.input)
    
    while(1):
        if "mdtrj" in arg.args.input:
            arg.args.mode = 'q'
            break
        elif "xsf" in arg.args.input:
            arg.args.mode = 'c'
            break
        elif "gro" in arg.args.input:
            arg.args.mode = 'g'
            break
        elif "lammpstrj" in arg.args.input:
            arg.args.mode = 'l'
            break
        else:
            kind = input("Iutput format. CONQUEST (c), Gromacs (g), LAMMPS (l), or QE (q): ")
            if kind.lower() == 'q':
                arg.input_type = "mdtrj"
            elif kind.lower() == 'c':
                arg.input_type = "xsf"
            elif kind.lower() == 'g':
                arg.input_type = "gro"
            elif kind.lower() == 'l':
                arg.input_type = "lammpstrj"
            else:
                print("Warning: Input 'c', 'g', 'l', or 'q'.")
                continue
            arg.args.mode = kind.lower()
            arg.input_check()
            break
    if arg.args.mode == 'q':
        if arg.args.elements is None:
            arg.args.elements = arg.input_file_check(getcwd(), "QE.in")
        qe = QE(arg)
        qe.read_all()
        xdatcar = XDATCAR(qe)
        xdatcar.write_all()
    elif arg.args.mode == 'c':
        conquest = Conquest(arg)
        conquest.read_all()
        xdatcar = XDATCAR(conquest)
        xdatcar.write_all()
    elif arg.args.mode == 'g':
        atom_step = AtomStep_Trj(arg)
        with open(arg.args.input) as read_file:
            with open(arg.args.output, "w") as write_file:
                molecule_count = -1
                step = 1
                first_flag = True
                for line in read_file:
                    if "t" in line and "step" in line:
                        title = line.split()[0]
                        line = next(read_file)
                        atom_sum = int(line.split()[0])
                        atom_info = defaultdict(list)
                        elements_mole_dict = defaultdict(list)
                        if first_flag:
                            former_mole_kind = ""
                        for atom in range(atom_sum):
                            line = next(read_file)
                            element = compile(r"([A-Za-z]+)(\d+)*").match(line.split()[1]).groups()[0]
                            atom_info[atom] = [element, float(line.split()[3])*10, float(line.split()[4])*10, float(line.split()[5])*10]
                            if first_flag:
                                mole_kind = compile(r"(\d+)([A-Za-z0-9]+)").match(line.split()[0]).groups()
                                if f"{mole_kind[1]}_{mole_kind[0]}" != former_mole_kind:
                                    former_mole_kind = f"{mole_kind[1]}_{mole_kind[0]}"
                                    atom_step.molecules.molecule_kind.append(f"{mole_kind[1]}_{mole_kind[0]}")
                                elements_mole_dict[element].append(f"{mole_kind[1]}_{mole_kind[0]}")
                        line = next(read_file)
                        lat_mat = zeros((3, 3))
                        if len(line.split()) == 3:
                            lat_mat[0, 0] = float(line.split()[0])*10; lat_mat[1, 1] = float(line.split()[1])*10; lat_mat[2, 2] = float(line.split()[2])*10
                        elif len(line.split()) == 9:
                            lat_mat[0, 0] = float(line.split()[0])*10; lat_mat[1, 1] = float(line.split()[1])*10; lat_mat[2, 2] = float(line.split()[2])*10
                            lat_mat[0, 1] = float(line.split()[3])*10; lat_mat[0, 2] = float(line.split()[4])*10
                            lat_mat[1, 0] = float(line.split()[5])*10; lat_mat[1, 2] = float(line.split()[6])*10
                            lat_mat[2, 0] = float(line.split()[7])*10; lat_mat[2, 1] = float(line.split()[8])*10
                        atom_step.lattice = lat_mat
                        ele_num, index_dict = rearange(atom_info)
                        if first_flag or atom_step.Lattice.NpT_flag:
                            write_file.write(f"{title}\n 1.0\n")
                            lat = atom_step.lattice if atom_step.lattice.ndim == 2 else atom_step.lattice[-1]
                            for arr in lat:
                                write_file.write(f"    {arr[0]:.10f} {arr[1]:.10f} {arr[2]:.10f}\n")
                            write_file.write("".join(f"   {idx}" for idx in ele_num.keys()) + '\n')
                            write_file.write("".join(f"   {idx}" for idx in ele_num.values()) + '\n')
                        write_file.write(f"Direct configuration=    {step}\n")
                        positions = []
                        for i in range(atom_sum):
                            positions.append(array(atom_info[index_dict[i]][1:]))
                        positions = atom_step.c2f(array(positions), lat)
                        for i in range(atom_sum):
                            write_file.write(f"   {positions[i][0]:.8f}    {positions[i][1]:.8f}    {positions[i][2]:.8f}\n")
                        if first_flag:
                            first_flag = False
                            for molecule_count, mole_kind in enumerate(atom_step.molecules.molecule_kind):
                                counter = 0
                                for values in elements_mole_dict.values():
                                    for value in values:
                                        if value == mole_kind:
                                            atom_step.molecules.molecule_dictionary[molecule_count].append(counter)
                                        counter += 1
                            if not path.isfile("molecules.csv"):
                                with open("molecules.csv", "w") as write_molecule_file:
                                    former_mole_kind = ""
                                    for mol_idx, mole_kind in enumerate(atom_step.molecules.molecule_kind):
                                        sp = mole_kind.split("_")
                                        if sp[0] != former_mole_kind:
                                            former_mole_kind = sp[0]
                                            write_molecule_file.write(f"mol {sp[0]}\n")
                                        write_molecule_file.write("".join([str(atom) if idx == 0 else f",{atom}" for idx, atom in enumerate(atom_step.molecules.molecule_dictionary[mol_idx])]) + "\n")
                        step += 1
    elif arg.args.mode == 'l':
        atom_step = AtomStep_Trj(arg)
        with open(arg.args.input) as read_file:
            first_flag = True
            atom_info = {}
            step = 1
            with open(arg.args.output, "w") as write_file:
                for line in read_file:
                    if "ITEM: TIMESTEP" in line:
                        line = next(read_file)
                        line = next(read_file)
                        line = next(read_file)
                        atom_sum = int(line.split()[0])
                        line = next(read_file)
                        sp = next(read_file).split()
                        tmp = zeros((3, 3))
                        i = 0
                        cell_lower_X, cell_upper_X = float(sp[0]), float(sp[1])
                        tmp[i, i] = cell_upper_X - cell_lower_X
                        sp = next(read_file).split()
                        i += 1
                        cell_lower_Y, cell_upper_Y = float(sp[0]), float(sp[1])
                        tmp[i, i] = cell_upper_Y - cell_lower_Y
                        sp = next(read_file).split()
                        i += 1
                        cell_lower_Z, cell_upper_Z = float(sp[0]), float(sp[1])
                        tmp[i, i] = cell_upper_Z - cell_lower_Z
                        atom_step.lattice = tmp
                        elements_mole_dict = defaultdict(list)
                        line = next(read_file)
                        for i in range(atom_sum):
                            sp = next(read_file).split()
                            atom_info[int(sp[0])] = [sp[3], float(sp[5])+cell_lower_X, float(sp[6])+cell_lower_Y, float(sp[7])+cell_lower_Z]
                            elements_mole_dict[int(sp[1])].append(int(sp[0]))
                        if first_flag or atom_step.Lattice.NpT_flag:
                            write_file.write("trajectory\n  1.0\n")
                            lat = atom_step.lattice if atom_step.lattice.ndim == 2 else atom_step.lattice[-1]
                            for l in lat:
                                write_file.write(f"   {l[0]:.10f}   {l[1]:.10f}   {l[2]:.10f}\n")
                            if first_flag:
                                first_flag = False
                                ele_num, index_dict = rearange(atom_info)
                            write_file.write("".join(f"   {idx}" for idx in ele_num.keys()) + '\n')
                            write_file.write("".join(f"   {idx}" for idx in ele_num.values()) + '\n')
                            if not path.isfile("molecules.csv"):
                                with open("molecules.csv", "w") as write_molecule_file:
                                    index_dict_new = {val: key for key, val in index_dict.items()}
                                    for mol_idx in sorted(elements_mole_dict.keys()):
                                        mole_atoms = [index_dict_new[atom] if idx == 0 else index_dict_new[atom] for idx, atom in enumerate(elements_mole_dict[mol_idx])]
                                        write_molecule_file.write("".join([str(atom) if idx == 0 else f",{atom}" for idx, atom in enumerate(sorted(mole_atoms))]) + "\n")
                        write_file.write(f"Direct configuration=    {step}\n")
                        positions = []
                        for i in range(atom_sum):
                            positions.append(array(atom_info[index_dict[i]][1:]))
                        positions = atom_step.c2f(array(positions), lat)
                        for i in range(atom_sum):
                            write_file.write(f"   {positions[i][0]:.8f}    {positions[i][1]:.8f}    {positions[i][2]:.8f}\n")
                        step += 1
    print("Done!")
