from re import compile
from numpy import zeros, array
import numpy as np
from os import path
from kit.fundamental import *
from kit.accelerate import *

class QE(Position, Periodic, Trajectory):
    def __init__(self, args=None, steps=None, atoms=None, **kwargs):
        super().__init__(input_obj=args, steps=steps, atoms=atoms, **kwargs)
        if hasattr(self, "_args") and hasattr(self._args, "input") and "mdtrj" in self._args.input:
            self.read_elements()
    @property
    def bridge(self):
        pass
    @bridge.setter
    def bridge(self, software):
        self.args, self.atom_step, self.elements = software.args, software.atom_step, software.elements
    def read_elements(self):
        assert path.isfile(self._args.elements), "Provide information about elements. They could be in 'QE.in'."
        with open(self._args.elements) as read_file:
            ele_flag = False
            for line in read_file:
                if "ATOMIC_POSITIONS" in line:
                    ele_flag = True
                elif ele_flag and len(line.split()) < 4:
                    ele_flag = False
                elif ele_flag:
                    self._elements.append(line.split()[0])
            self.__atom_sum = len(self._elements)
    def read_all(self, atoms="all", steps="all"):
        if atoms == "all":
            self.atoms.bridge = self
            self.atoms.put(atoms)
        else:
            self.atoms = atoms
        atoms = (atoms if atoms == "all" else self._AtomStep.atom_list)
        if atoms != "all":
            atoms_dict = {atom: atom in atoms for atom in range(self.__atom_sum)}
        if steps == "all":
            self.steps.bridge = self
            self.steps.put(steps)
        else:
            self.steps = steps
            self.cutoff_step = max(steps.get(slice_flatten=True))
            if self._cutoff_step > 1000:
                print("Completed steps:")
        steps = (steps if steps == "all" else self._AtomStep.steps.get(slice_flatten=True))
        if steps != "all":
            steps_dict = {step: step in steps for step in range(1, self.steps.total_steps+1)}
        with open(self._args.input) as read_file:
            self._AtomStep.cartesian_flag, self._AtomStep.fractional_flag = True, False
            
            for line in read_file:
                self._step_counter += 1
                if not (steps == "all" or steps_dict[self._step_counter]):
                    for _ in range(self.__atom_sum+5):
                        line = next(read_file)
                else:
                    if self._step_counter == 1:
                        line = next(read_file)
                    lat_mat = zeros((3, 3))
                    for i in range(3):
                        lat_mat[i] = [float(_) for _ in line.split()]
                        line = next(read_file)
                    self._AtomStep.lattice = lat_mat

                    for atom in range(self.__atom_sum):
                        if (atoms == "all" or atoms_dict[atom]):
                            self._AtomStep.put([float(_) for _ in line.split()], atom)
                        line = next(read_file)

                    line = next(read_file, None)
                    if line is None:
                        break
                    while(len(line.split()) != 3):
                        line = next(read_file)

                    if self._step_counter > self._cutoff_step and self._cutoff_step != -1:
                        break
            
            self._AtomStep.lock()
    def write_all(self):
        self._args.output = self._title
        if "mdtrj" not in self._args.output:
            self._args.output += ".mdtrj"
        with open(self._args.output, "w") as write_file:
            for idx, cart_pos in enumerate(self._AtomStep.cartesian_position):
                write_file.write("{:.12f}\t{:.12f}\t{:.12f}\n".write(0, 0, 0))
                if self._AtomStep.lattice.ndim == 3:
                    for _ in self._AtomStep.lattice[idx]:
                        write_file.write(f"{_[0]:.12f}\t{_[1]:.12f}\t{_[2]:.12f}\n")
                else:
                    for _ in self._AtomStep.lattice:
                        write_file.write(f"{_[0]:.12f}\t{_[1]:.12f}\t{_[2]:.12f}\n")
                for position in cart_pos:
                    write_file.write(f"{position[0]:.12f}\t{position[1]:.12f}\t{position[2]:.12f}\n")
                write_file.write('\n')

class QE_Single_Point(Single_Point, Fundamental):
    def __init__(self, args=None, atoms=None, **kwargs):
        super().__init__(input_obj=args, atoms=atoms, **kwargs)
        if hasattr(self, "_args") and hasattr(self._args, "elements"):
            self.read_elements()
    @property
    def bridge(self):
        pass
    @bridge.setter
    def bridge(self, software):
        self.args, self.elements, self.atom_step, self.relax, self.relax_border = software.args, software.elements, software.atom_step, software.relax, software.relax_border
    def read_elements(self):
        with open(self._args.input) as read_file:
            for line in read_file:
                if "ATOMIC_POSITIONS" in line.split():
                    while(1):
                        line = next(read_file)
                        if len(line.split()) < 4:
                            break
                        self._elements.append(line.split()[0])
                    self._elements = tuple(self._elements)
        return self._elements
    def read_all(self, atoms="all"):
        if atoms == "all":
            self.atoms.bridge = self
            self.atoms.put(atoms)
        else:
            self.atoms = atoms
        atoms = (atoms if atoms == "all" else self._AtomStep.atom_list)
        if atoms != "all":
            atoms_dict = {atom: atom in atoms for atom in range(len(self._elements))}
        with open(self._args.input) as read_file:
            for line in read_file:
                if "ATOMIC_POSITIONS" in line.split():
                    if "angstrom" in line.split() or "bohr" in line.split():
                        self._AtomStep.fractional_flag, self._AtomStep.cartesian_flag = False, True
                    elif "crystal" in line.split():
                        self._AtomStep.fractional_flag, self._AtomStep.cartesian_flag = True, False
                    else:
                        raise Exception("Only access 'bohr', 'angstrom', or 'crystal' in 'ATOMIC_POSITIONS' row.")
                    atom = 0
                    while(1):
                        line = next(read_file)
                        if len(line.split()) < 4:
                            break
                        sp = line.split()
                        if (atoms == "all" or atoms_dict[atom]):
                            self._AtomStep.put([float(_) for _ in sp[1:4]], atom)
                        atom += 1
                        if len(sp) > 4:
                            self._relax_flag = True
                            relax_direction = ""
                            for _ in range(4, 7):
                                relax_direction += ("T" if sp[_] == "1" else "F")
                            self._relax_border.append(relax_direction)
                        else:
                            self._relax_border.append("TTT")
                elif "CELL_PARAMETERS" in line.split():
                    if "bohr" in line.split():
                        self._scaling = 0.529177249
                    lat_mat = zeros((3, 3))
                    for _ in range(3):
                        line = next(read_file)
                        lat_mat[_] = [float(_) for _ in line.split()]
                    self._AtomStep.lattice = lat_mat
                elif "prefix" in line.split():
                    self._title = line.replace("'", "").replace(",", "").replace("=", "").split()[1]
            self._AtomStep.lock()
    def write_all(self, unit="crystal"):
        self._AtomStep.lattice = self._scaling*self._AtomStep.lattice
        if unit == "angstrom" or unit == "bohr":
            position = self._AtomStep.cartesian_position = self._scaling*self._AtomStep.cartesian_position
        elif unit == "crystal":
            position = self._AtomStep.fractional_position
        with open(self._args.output, "w") as write_file:
            write_file.write(f"ATOMIC_POSITIONS {unit}\n")
            for idx, element in enumerate(self._AtomStep.atoms.elements):
                write_file.write(f"{element} {position[idx][0]:.8f} {position[idx][1]:.8f} {position[idx][2]:.8f}")
                if self._relax_flag and isinstance(self._relax_border, list):
                    for _ in range(len(self._relax_border[idx])):
                        write_file.write(f" {'1' if self._relax_border[idx][_] == 'T' else '0'}")
                elif self._relax_flag and isinstance(self._relax_border, float):
                    write_file.write((" 0 0 0" if position[idx][2] <= self._relax_border else " 1 1 1"))
                write_file.write('\n')
            write_file.write("\n\n")
            write_file.write(f"CELL_PARAMETERS {unit}\n")
            for _ in self._AtomStep.lattice:
                write_file.write(f"{_[0]:.10f}\t{_[1]:.10f}\t{_[2]:.10f}\n")

class Conquest(Position, Periodic, Trajectory):
    def __init__(self, args=None, steps=None, atoms=None, **kwargs):
        super().__init__(input_obj=args, steps=steps, atoms=atoms, **kwargs)
        if hasattr(self, "_args") and hasattr(self._args, "elements"):
            self.read_elements()
        if hasattr(args, "bridge"):
            self.bridge = args
    @property
    def bridge(self):
        pass
    @bridge.setter
    def bridge(self, software):
        self.args, self.atom_step, self.elements = software.args, software.atom_step, software.elements
    def read_elements(self):
        from os import popen
        with popen(f"head {self._args.input}") as command:
            command = iter(command)
            for line in command:
                if "PRIMCOORD" in line.upper():
                    line = next(command)
                    self.__atom_sum = int(line.split()[0])
                    break
        with popen(f"head -{self.__atom_sum+10} {self._args.input}") as command:
            command = iter(command)
            for line in command:
                if "PRIMCOORD" in line.upper():
                    line = next(command)
                    line = next(command)
                    for _ in range(self.__atom_sum):
                        self._elements.append(line.split()[0])
                        line = next(command)
                    self._elements = tuple(self._elements)
                    break
    def read_all(self, steps="all", atoms="all"):
        if atoms == "all":
            self.atoms.bridge = self
            self.atoms.put(atoms)
        else:
            self.atoms = atoms
        atoms = (atoms if atoms == "all" else self._AtomStep.atom_list)
        if atoms != "all":
            atoms_dict = {atom: atom in atoms for atom in range(len(self._elements))}
        if steps == "all":
            self.steps.bridge = self
            self.steps.put(steps)
        else:
            self.steps = steps
            self.cutoff_step = max(steps.get(slice_flatten=True))
            if self._cutoff_step > 1000:
                print("Completed steps:")
        steps = (steps if steps == "all" else self._AtomStep.steps.get(slice_flatten=True))
        if steps != "all":
            steps_dict = {step: step in steps for step in range(1, self.steps.total_steps+1)}
        with open(self._args.input) as read_file:
            for line in read_file:
                if "PRIMVEC" in line.upper():
                    self._step_counter += 1

                    lat_mat = zeros((3, 3))
                    for i in range(3):
                        line = next(read_file)
                        lat_mat[i] = [float(_) for _ in line.split()]
                    self._AtomStep.lattice = lat_mat

                    line = next(read_file)
                    line = next(read_file)
                    if atoms != "all":
                        atoms_dict = {atom: atom in atoms for atom in range(self.__atom_sum)}
                    self._AtomStep.cartesian_flag, self._AtomStep.fractional_flag = True, False

                    if not (steps == "all" or steps_dict[self._step_counter]):
                        for _ in range(self.__atom_sum):
                            line = next(read_file)
                    elif len(self._elements) == 0:
                        for atom in range(self.__atom_sum):
                            line = next(read_file)
                            self._elements.append(line.split()[0])
                            if (atoms == "all" or atoms_dict[atom]):
                                self._AtomStep.put([float(_) for _ in line.split()[1:]], atom)
                        self._elements = tuple(self._elements)
                    else:
                        for atom in range(self.__atom_sum):
                            line = next(read_file)
                            if (atoms == "all" or atoms_dict[atom]):
                                self._AtomStep.put([float(_) for _ in line.split()[1:]], atom)
                
                if self._step_counter > self._cutoff_step and self._cutoff_step != -1:
                    break
            self._AtomStep.lock()
    def write_all(self):
        with open(self._args.output, "w") as write_file:
            for idx, cart_pos in enumerate(self._AtomStep.cartesian_position):
                write_file.write("CRYSTAL\n")
                write_file.write(f"PRIMVEC     \t{idx}\n")
                if self._AtomStep.lattice.ndim == 3:
                    for lat_arr in self._AtomStep.lattice[idx]:
                        write_file.write(f"{lat_arr[0]:.12f}\t{lat_arr[1]:.12f}\t{lat_arr[2]:.12f}\n")
                else:
                    for lat_arr in self._AtomStep.lattice:
                        write_file.write(f"{lat_arr[0]:.12f}\t{lat_arr[1]:.12f}\t{lat_arr[2]:.12f}\n")
                write_file.write(f"PRIMCOORD   \t{idx}\n")
                write_file.write(f"    {len(self._AtomStep.elements)}  1\n")
                for _, position in enumerate(cart_pos):
                    write_file.write(f"  {self._AtomStep.elements[_]}\t{position[0]:.10f}\t{position[1]:.10f}\t{position[2]:.10f}\n")
                write_file.write('\n')

class Conquest_Single_Point(Single_Point, Fundamental):
    def __init__(self, args=None, atoms=None, **kwargs):
        super().__init__(input_obj=args, atoms=atoms, **kwargs)
        if hasattr(self, "_args") and hasattr(self._args, "elements"):
            self.read_elements()
        if hasattr(args, "bridge"):
            self.bridge = args
    @property
    def bridge(self):
        pass
    @bridge.setter
    def bridge(self, software):
        self.args, self.elements, self.atom_step, self.relax, self.relax_border = software.args, software.elements, software.atom_step, software.relax, software.relax_border
    def read_elements(self):
        with open(self._args.elements) as read_file:
            element_order, elements = {}, {}
            for line in read_file:
                if "ChemicalSpeciesLabel" in line:
                    while(1):
                        line = next(read_file)
                        if "endblock" in line:
                            break
                        sp = line.split()
                        element_order[sp[0]] = sp[2]
                        elements[sp[2]] = 0
        with open(self._args.input) as read_file:
            row = 1
            for line in read_file:
                if row > 4:
                    sp = line.split()
                    elements[element_order[sp[3]]] += 1
                row += 1
        ele = [element_order[str(i)] for i in range(1, len(element_order.keys())+1)]
        num = [elements[element_order[str(i)]] for i in range(1, len(element_order.keys())+1)]
        self._AtomStep.elements = Convert.eleNum2elements(ele, num)
        self._elements = tuple(self._AtomStep.elements)
        self.__atom_sum = len(self._elements)
        return self._AtomStep.elements
    def read_all(self, atoms="all"):
        if atoms == "all":
            self.atoms.bridge = self
            self.atoms.put(atoms)
        else:
            self.atoms = atoms
        atoms = (atoms if atoms == "all" else self._AtomStep.atom_list)
        if atoms != "all":
            atoms_dict = {atom: atom in atoms for atom in range(len(self._elements))}
        with open(self._args.elements) as read_file:
            element_order, elements = {}, {}
            for line in read_file:
                if "ChemicalSpeciesLabel" in line:
                    while(1):
                        line = next(read_file)
                        if "endblock" in line:
                            break
                        sp = line.split()
                        element_order[sp[0]] = sp[2]
                        elements[sp[2]] = 0
                    break
                elif "IO.Title" in line.split():
                    if len(line.split()) == 1:
                        self._title = self._args.input.split('.')[0]
                    else:
                        self._title = line.split()[1]
        with open(self._args.input) as read_file:
            self._relax_border = []
            for line in read_file:
                lattice_mat = zeros((3, 3))
                for i in range(3):
                    lattice_mat[i] = [float(_) for _ in line.split()]
                    line = next(read_file)
                self._AtomStep.lattice = lattice_mat

                self._AtomStep.fractional_flag = True

                for atom in range(self.__atom_sum):
                    sp = next(read_file).split()
                    if not (atoms == "all" or atoms_dict[atom]):
                        continue
                    elements[element_order[sp[3]]] += 1
                    self._AtomStep.put([float(_) for _ in sp[:3]], atom)
                    if len(sp) > 4:
                        self._relax_flag = True
                        self._relax_border.append(f"{sp[4]}{sp[5]}{sp[6]}")
        self._AtomStep.lock()
        ele = [element_order[str(_)] for _ in range(1, len(element_order.keys())+1)]
        num = [elements[element_order[str(_)]] for _ in range(1, len(element_order.keys())+1)]
        self._elements = tuple(Convert.eleNum2elements(ele, num))
        if atoms != "all":
            self._AtomStep.atoms.elements = [self._AtomStep.atoms.elements[atom] for atom in atoms]
    def write_all(self):
        with open(self._args.output, "w") as write_file:
            for lat_arr in self._scaling*self._AtomStep.lattice:
                write_file.write(f"{lat_arr[0]:.10f}\t{lat_arr[1]:.10f}\t{lat_arr[2]:.10f}\n")
            write_file.write(f"{len(self._AtomStep.elements)}\n")
            element_num = 1
            for idx, position in enumerate(self._AtomStep.fractional_position):
                if idx > 0 and self._AtomStep.elements[idx-1] != self._AtomStep.elements[idx]:
                    element_num += 1
                write_file.write(f"{position[0]:.8f} {position[1]:.8f} {position[2]:.8f} {element_num} ")
                if self._relax_flag and isinstance(self._relax_border, list):
                    write_file.write(f" {self._relax_border[idx][0]} {self._relax_border[idx][1]} {self._relax_border[idx][2]}\n")
                elif self._relax_flag and isinstance(self._relax_border, float):
                    write_file.write((" F F F\n" if position[idx][2] <= self._relax_border else " T T T\n"))
                elif not self._relax_flag:
                    write_file.write('T T T\n')

class Gromacs(Position, Periodic, Trajectory):
    def __init__(self, args=None, steps=None, atoms=None, **kwargs):
        super().__init__(input_obj=args, steps=steps, atoms=atoms, **kwargs)
        if hasattr(self, "_args") and hasattr(self._args, "input") and "gro" in self._args.input:
            self.read_elements()
        if hasattr(args, "bridge"):
            self.bridge = args
    @property
    def bridge(self):
        pass
    @bridge.setter
    def bridge(self, software):
        self.args, self.elements, self.atom_step = software.args, software.elements, software.atom_step
    def read_elements(self):
        from os import popen
        with popen(f"sed -n '2p'  {self._args.input}") as command:
            self.__atom_sum = int(command.read())
        with popen(f"head -{self.__atom_sum+2} {self._args.input}") as command:
            for line in command:
                line = next(iter(command))
                for _ in range(self.__atom_sum):
                    line = next(iter(command))
                    self._elements.append(compile(r"([A-Za-z]+)").match(line.split()[1]).groups()[0])
            self._AtomStep.molecules.elements = self._AtomStep.atoms.elements = self._elements
        self._elements = tuple(self._elements)
        return self._elements
    def read_all(self, atoms="all", steps="all", rearange=True):
        if atoms == "all":
            self.atoms.bridge = self
            self.atoms.put("all")
        else:
            self.atoms = atoms
        atoms = (atoms if atoms == "all" else self._AtomStep.atom_list)
        if atoms != "all":
            atoms_dict = {atom: atom in atoms for atom in range(len(self._elements))}
        if steps == "all":
            self.steps.bridge = self
            self.steps.put("all")
        else:
            self.steps = steps
            self.cutoff_step = max(steps.get())
        steps = (steps if steps == "all" else self._AtomStep.steps.get())
        if steps != "all":
            steps_dict = {step: step in steps for step in range(1, self.steps.total_steps+1)}
        with open(self._args.input) as read_file:
            if rearange:
                cart_pos_steps = []
            self._AtomStep.cartesian_flag, self._AtomStep.fractional_flag = True, False
            molecule_count = -1
            for line in read_file:
                if "t" in line and "step" in line:
                    self._step_counter += 1
                    line = next(read_file)
                    if rearange:
                        elements_dict, elements_mole_dict = defaultdict(list), defaultdict(list)
                    if not (steps == "all" or steps_dict[self._step_counter]):
                        for _ in range(self.__atom_sum):
                            line = next(read_file)
                    if atoms != "all":
                        atoms_dict = {atom: atom in atoms for atom in range(self.__atom_sum)}
                    else:
                        if len(self._AtomStep.molecules.molecule_kind) == 0:
                            former_mole_kind = ""
                        for atom in range(self.__atom_sum):
                            line = next(read_file)
                            element = compile(r"([A-Za-z]+)(\d+)*").match(line.split()[1]).groups()[0]
                            if not (atoms == "all" or atoms_dict[atom]):
                                continue
                            if self._step_counter == 1:
                                mole_kind = compile(r"(\d+)([A-Za-z0-9]+)(\d+)?").match(line.split()[0]).groups()
                                if f"{mole_kind[1]}_{mole_kind[0]}" != former_mole_kind:
                                    former_mole_kind = f"{mole_kind[1]}_{mole_kind[0]}"
                                    self._AtomStep.molecules.molecule_kind.append(f"{mole_kind[1]}_{mole_kind[0]}")
                                    if rearange:
                                        elements_mole_dict[element].append(f"{mole_kind[1]}_{mole_kind[0]}")
                                    else:
                                        molecule_count += 1
                                        self._AtomStep.molecules.molecule_dictionary[molecule_count].append(atom)
                            if rearange:
                                tmp = [atom]
                                tmp.extend([float(_)*10 for _ in line.split()[3:6]])
                                elements_dict[element].append(tmp)
                            else:
                                self._AtomStep.atoms.elements.append(element)
                                self._AtomStep.put([float(_)*10 for _ in line.split()[3:6]], atom)
                    line = next(read_file)
                    lat_mat = zeros((3, 3))
                    if len(line.split()) == 3:
                        lat_mat[0][0] = float(line.split()[0])*10; lat_mat[1][1] = float(line.split()[1])*10; lat_mat[2][2] = float(line.split()[2])*10
                    elif len(line.split()) == 9:
                        lat_mat[0][0] = float(line.split()[0])*10; lat_mat[1][1] = float(line.split()[1])*10; lat_mat[2][2] = float(line.split()[2])*10
                        lat_mat[0][1] = float(line.split()[3])*10; lat_mat[0][2] = float(line.split()[4])*10
                        lat_mat[1][0] = float(line.split()[5])*10; lat_mat[1][2] = float(line.split()[6])*10
                        lat_mat[2][0] = float(line.split()[7])*10; lat_mat[2][1] = float(line.split()[8])*10
                    self._AtomStep.lattice = lat_mat
                    if rearange:
                        cart_pos = []
                        if self._step_counter == 1:
                            self._elements, elements = [], []
                        num = 0
                        for element, values in elements_dict.items():
                            for value in values:
                                if self._step_counter == 1:
                                    self._elements.append(element)
                                    elements.append(self._AtomStep.atoms.elements[value[0]])
                                cart_pos.append(value[1:])
                                num += 1
                        if self._step_counter == 1:
                            self._elements = tuple(self._elements)
                            self._AtomStep.elements = elements
                        
                        if len(self._AtomStep.molecules.molecule_dictionary) == 0:
                            for molecule_count, mole_kind in enumerate(self._AtomStep.molecules.molecule_kind):
                                counter = 0
                                for values in elements_mole_dict.values():
                                    for value in values:
                                        if value == mole_kind:
                                            self._AtomStep.molecules.molecule_dictionary[molecule_count].append(counter)
                                        counter += 1

                        cart_pos_steps.append(cart_pos)
            if not path.isfile("molecules.csv"):
                with open("molecules.csv", "w") as write_file:
                    former_mole_kind = ""
                    for mol_idx, mole_kind in enumerate(self._AtomStep.molecules.molecule_kind):
                        sp = mole_kind.split("_")
                        if sp[0] != former_mole_kind:
                            former_mole_kind = sp[0]
                            write_file.write(f"mol {sp[0]}\n")
                        write_file.write("".join([str(atom) if idx == 0 else f",{atom}" for idx, atom in enumerate(self._AtomStep.molecules.molecule_dictionary[mol_idx])]) + "\n")
            if rearange:
                self._AtomStep.cartesian_position = array(cart_pos_steps)
            else:
                self._AtomStep.lock()

class Gromacs_Single_Point(Single_Point):
    def __init__(self, args=None, atoms=None, **kwargs):
        super().__init__(input_obj=args, atoms=atoms, **kwargs)
        if hasattr(self, "_args") and hasattr(self._args, "input") and "gro" in self._args.input:
            if self._args.input is not None and path.isfile(self._args.input):
                self.read_elements()
        if hasattr(args, "bridge"):
            self.bridge = args
    @property
    def bridge(self):
        pass
    @bridge.setter
    def bridge(self, software):
        self.args, self.elements, self.atom_step = software.args, software.elements, software.atom_step
    def read_elements(self):
        with open(self._args.input) as read_file:
            for line in read_file:
                line = next(read_file)
                self.__atom_sum = int(line.split()[0])
                for _ in range(self.__atom_sum):
                    sp = next(read_file).split()
                    if len(sp) == 5 and compile(r"[A-Za-z]+\d+").match(sp[1]) is not None:
                        tmp = compile(r"([A-Za-z]+)(\d+)").match(sp[1]).groups()
                        sp[1] = tmp[0]; sp.insert(2, tmp[1])
                    self._elements.append(sp[1])
                break
            self._elements = tuple(self._elements)
        return self._elements
    def read_all(self, atoms="all", rearange=True):
        if atoms == "all":
            self.atoms.bridge = self
            self.atoms.put(atoms)
        else:
            self.atoms = atoms
        atoms = (atoms if atoms == "all" else self._AtomStep.atom_list)
        if atoms != "all":
            atoms_dict = {atom: atom in atoms for atom in range(len(self._elements))}
        self._AtomStep.cartesian_flag = True
        with open(self._args.elements) as read_file:
            mole_kinds = []
            for line in read_file:
                if "molecules" in line:
                    while(1):
                        line = next(read_file)
                        if len(line.split()) == 2:
                            mole_kinds.append(line.split()[0])
                        elif len(mole_kinds) == 0:
                            pass
                        else:
                            break
        if rearange:
            elements_dict, elements_mole_dict = defaultdict(list), defaultdict(list)
        mol_idx, former_mol_idx, mole_kinds_idx = 0, 1000000, -1
        with open(self._args.input) as read_file:
            for line in read_file:
                line = next(read_file)
                for atom in range(self.__atom_sum):
                    sp = next(read_file).split()
                    if not (atoms == "all" or atoms_dict[atom]):
                        continue
                    if len(sp) == 5 and compile(r"[A-Za-z]+\d+").match(sp[1]) is not None:
                        tmp = compile(r"([A-Za-z]+)(\d+)").match(sp[1]).groups()
                        sp[1] = tmp[0]; sp.insert(2, tmp[1])
                    if former_mol_idx > int(sp[0]):
                        former_mol_idx = int(sp[0])
                        mole_kinds_idx += 1
                        mol_idx += 1
                    elif former_mol_idx != int(sp[0]):
                        former_mol_idx = int(sp[0])
                        mol_idx += 1
                    mole_kind = f"{mole_kinds[mole_kinds_idx]}_{sp[0]}"
                    if mole_kind not in self._AtomStep.molecules.molecule_kind:
                        self._AtomStep.molecules.molecule_kind.append(mole_kind)
                    if rearange:
                        tmp = [atom]
                        tmp.extend([float(_) for _ in sp[3:6]])
                        elements_dict[sp[1]].append(tmp)
                        elements_mole_dict[sp[1]].append(mole_kind)
                    else:
                        self._AtomStep.put([float(_) for _ in sp[3:6]], atom)
                        self._AtomStep.molecules.molecule_dictionary[mol_idx].append(atom)
                if rearange:
                    self._elements, elements, cart_pos = [], [], []
                    num = 0
                    for element, atom_positions in elements_dict.items():
                        for idx, atom_position in enumerate(atom_positions):
                            cart_pos.append(atom_position[1:])
                            self._elements.append(element)
                            elements.append(self._AtomStep.atoms.elements[atom_position[0]])
                            mol_idx = self._AtomStep.molecules.molecule_kind.index(elements_mole_dict[element][idx])
                            self._AtomStep.molecules.molecule_dictionary[mol_idx].append(num)
                            num += 1
                    self._AtomStep.elements, self._AtomStep.cartesian_position = elements, array(cart_pos)
                    self._elements = tuple(self._elements)
                else:
                    self._AtomStep.lock()
                line = next(read_file)
                lat_mat = zeros((3, 3))
                if len(line.split()) == 3:
                    lat_mat[0][0] = float(line.split()[0]); lat_mat[1][1] = float(line.split()[1]); lat_mat[2][2] = float(line.split()[2])
                elif len(line.split()) == 9:
                    lat_mat[0][0] = float(line.split()[0]); lat_mat[1][1] = float(line.split()[1]); lat_mat[2][2] = float(line.split()[2])
                    lat_mat[0][1] = float(line.split()[3]); lat_mat[0][2] = float(line.split()[4])
                    lat_mat[1][0] = float(line.split()[5]); lat_mat[1][2] = float(line.split()[6])
                    lat_mat[2][0] = float(line.split()[7]); lat_mat[2][1] = float(line.split()[8])
                self._AtomStep.lattice = lat_mat
                if not path.isfile("molecules.csv"):
                    with open("molecules.csv", "w") as write_file:
                        former_mole_kind = ""
                        for mol_id, mole_kind in enumerate(self._AtomStep.molecules.molecule_kind):
                            if former_mole_kind != mole_kind.split("_")[0]:
                                write_file.write(f"mol {mole_kind.split('_')[0]}\n")
                                former_mole_kind = mole_kind.split("_")[0]
                            write_file.write("".join([str(atom) if idx == 0 else f",{atom}" for idx, atom in enumerate(self._AtomStep.molecules.molecule_dictionary[mol_id])]) + '\n')
    def write_all(self):
        with open(self._args.output, "w") as write_file:
            write_file.write(self._title + '\n')
            write_file.write(str(len(self._elements)) + '\n')
            num = 1
            for mol_id, mol_atoms in self._AtomStep.molecules.molecule_dictionary.items():
                for atom in mol_atoms:
                    write_file.write(f"{self._AtomStep.molecules.molecule_kind[mol_id].split('_')[1]:>5s}{self._elements[atom]:>10s}{num:5d}")
                    write_file.write(f"{self._AtomStep.cartesian_position[atom][0]*self._scaling:8.3f}{self._AtomStep.cartesian_position[atom][1]*self._scaling:8.3f}{self._AtomStep.cartesian_position[atom][2]*self._scaling:8.3f}\n")
                    num += 1
            write_file.write(f"  {self._AtomStep.lattice[0][0]*self._scaling:.5f}   {self._AtomStep.lattice[1][1]*self._scaling:.5f}   {self._AtomStep.lattice[2][2]*self._scaling:.5f}\n")

class Gaussian(Single_Point):
    def __init__(self, args=None, atoms=None, **kwargs):
        super().__init__(input_obj=args, atoms=atoms, **kwargs)
        self._relax_border = 0
        self.__charge = 0
        self.__multiplicity = 1
        self.__chk = None
        self.__method = " cam-b3lyp/6-311g geom=connectivity"
        if hasattr(self, "_args") and hasattr(self._args, "input") and "gjf" in self._args.input:
            self.read_elements()
        if hasattr(args, "bridge"):
            self.bridge = args
    @property
    def bridge(self):
        pass
    @property
    def fractional_position(self):
        return self._AtomStep.c2f(self._AtomStep.cartesian_position, self._AtomStep.lattice)
    @bridge.setter
    def bridge(self, software):
        self.args, self.title, self.elements, self.atom_step = software.args, software.title, software.elements, software.atom_step
    @fractional_position.setter
    def fractional_position(self, frac_pos):
        assert frac_pos.ndim == 2, "Can not provide trajectory data."
        self._AtomStep.fractional_position = array(frac_pos)
    def read_elements(self):
        with open(self._args.input) as read_file:
            for line in read_file:
                if compile(r"^-?\d+ \d+$").match(line.replace('\n', '')) is not None:
                    line = next(read_file)
                    while(line is not None and (len(line.split()) == 4)):
                        if line.split()[0] == "Tv":
                            line = next(read_file, None)
                        else:
                            self._elements.append(line.split()[0])
                            line = next(read_file, None)
                    break
        self._elements = tuple(self._elements)
        return self._elements
    def read_all(self, atoms="all", rearange=True):
        from collections import defaultdict
        if atoms == "all":
            self.atoms.bridge = self
            self.atoms.put(atoms)
        else:
            self.atoms = atoms
        atoms = (atoms if atoms == "all" else self._AtomStep.atom_list)
        if atoms != "all":
            atoms_dict = {atom: atom in atoms for atom in range(len(self._elements))}
        self._AtomStep.cartesian_flag = True
        with open(self._args.input) as read_file:
            num = 0
            lattice, cart_pos = [], []
            for line in read_file:
                if compile(r"^-?\d+ \d+$").match(line) is not None:
                    self.__charge = int(line.split()[0])
                    self.__multiplicity = int(line.split()[1])
                    line = next(read_file)
                    while(len(line.split()) != 4):
                        line = next(read_file)
                    if rearange:
                        elements_dict = defaultdict(list)
                    else:
                        elements = []
                    while(line is not None and (len(line.split()) == 4)):
                        if line.split()[0] == "Tv":
                            lattice.append([float(_) for _ in line.split()[1:]])
                        elif atoms == "all" or atoms_dict[num]:
                            if rearange:
                                elements_dict[line.split()[0]].append([float(_) for _ in line.split()[1:]])
                            else:
                                elements.append(line.split()[0])
                                self._AtomStep.put([float(_) for _ in line.split()[1:3]], num)
                        num += 1
                        line = next(read_file, None)
                    if lattice != []:
                        self._AtomStep.lattice = array(lattice)
                elif line[0] == '%':
                    self.__chk = line.split()[1:]
                elif line[0] == '#':
                    self.__method = line.split()[1:]
        
        if rearange:
            num = 0
            self._elements = []
            for key, values in elements_dict.items():
                for val in values:
                    self._elements.append(key)
                    cart_pos.append(val)
                    num += 1
            self._AtomStep.cartesian_position = array(cart_pos)
        else:
            self._AtomStep.lock()
        self._AtomStep.elements = self._elements
        self._elements = tuple(self._elements)
    def write_all(self):
        self._args.output = self._title
        if "gjf" not in self._title:
            self._args.output += ".gjf"
        with open(self._args.output, "w") as write_file:
            if self.__chk is not None:
                write_file.write(f"%{self.__chk}\n")
            write_file.write(f"#{self.__method}\n\n")
            write_file.write("Title Card Required\n\n")
            write_file.write(f"{self.__charge} {self.__multiplicity}\n")
            for element, cart_pos in zip(self._elements, self._AtomStep.cartesian_position):
                write_file.write(f"\t{element}\t{cart_pos[0]:.8f}\t{cart_pos[1]:.8f}\t{cart_pos[2]:.8f}\n")
            if list(self._AtomStep.lattice) != []:
                for lattice in self._AtomStep.lattice:
                    write_file.write(f"\tTv\t{lattice[0]:10f}\t{lattice[1]:.10f}\t{lattice[2]:.10f}\n")

class LAMMPS(Position, Periodic, Trajectory):
    def __init__(self, args=None, atoms=None, steps=None, **kwargs):
        super().__init__(input_obj=args, atoms=atoms, steps=steps, **kwargs)
        if hasattr(self, "_args") and hasattr(self._args, "elements"):
            self.read_elements()
        if hasattr(args, "bridge"):
            self.bridge = args
        self._AtomStep.steps = Step_LAMMPS()
    @property
    def bridge(self):
        pass
    @bridge.setter
    def bridge(self, software):
        self.args, self.title, self.elements, self.atom_step = software.args, software.title, software.elements, software.atom_step
    def read_elements(self):
        tmp = {}
        with open(self._args.input) as read_file:
            for _ in range(4):
                line = next(read_file)
            atom_sum = int(line.split()[0])
            for _ in range(5):
                line = next(read_file)
            for _ in range(atom_sum):
                line = next(read_file)
                sp = line.split()
                tmp[int(sp[0])-1] = sp[3]
            for key in sorted(tmp.keys()):
                self._elements.append(tmp[key])
        self._elements = tuple(self._elements)
        return self._elements
    def read_all(self, steps="all", atoms="all", rearange=True):
        from os import popen
        if atoms == "all":
            self.atoms.bridge = self
            self.atoms.put(atoms)
        else:
            self.atoms = atoms
        atoms = (atoms if atoms == "all" else self._AtomStep.atom_list)
        if atoms != "all":
            atoms_dict = {atom: atom in atoms for atom in range(len(self._elements))}
        if steps == "all":
            with popen(f"grep TIMESTEP {self._args.input} | wc") as command:
                self.steps.total_steps = int(command.read().split()[0])
            self.steps.put("all")
        else:
            self.steps = steps
            self.cutoff_step = max(steps.get())
        steps = (steps if steps == "all" else self._AtomStep.steps.get())
        if steps != "all":
            steps_dict = {step: step in steps for step in range(1, self.steps.total_steps+1)}
        self._AtomStep.cartesian_flag = True
        with open(self._args.input) as read_file:
            flag = True
            for line in read_file:
                if "TIMESTEP" in line:
                    line = next(read_file)
                    if self._step_counter > self._cutoff_step and self._cutoff_step != -1:
                        break
                    self._step_counter += 1
                    if not (steps == "all" or steps_dict[self._step_counter]):
                        line = next(read_file)
                        line = next(read_file)
                        for _ in range(int(line.split()[0])+5):
                            line = next(read_file)
                elif "BOX" in line:
                    lattice = zeros((3, 3))
                    for i in range(3):
                        sp = next(read_file).split()
                        lattice[i][i] = np.float32(sp[1]) - np.float32(sp[0])
                    self._AtomStep.lattice = lattice
                elif "NUMBER" in line:
                    sp = next(read_file).split()
                    if self._step_counter == 1:
                        atom_num = int(sp[0])
                elif "ATOMS" in line:
                    for atom in range(atom_num):
                        sp = next(read_file).split()
                        num = int(sp[0])-1
                        if not (atoms == "all" or atoms_dict[num]):
                            continue
                        self._AtomStep.put([np.float32(_) for _ in sp[-3:]], num)
                        if flag:
                            if atoms == "all":
                                self.atoms.put(atom)
                            self._AtomStep.atoms_info[num]["molecule"], self._AtomStep.atoms_info[num]["type"], self._AtomStep.atoms_info[num]["element"], self._AtomStep.atoms_info[num]["q"] = int(sp[1])-1, int(sp[2]), sp[3], float(sp[4])
                    if flag:
                        flag = False
            self._AtomStep.lock(dtype=np.float32)
    def write_all(self):
        if 'q' not in self._AtomStep.atoms_info[0].keys():
            for key in self._AtomStep.atoms_info.keys():
                self._AtomStep.atoms_info[key]["q"] = 0
        if 'type' not in self._AtomStep.atoms_info[0].keys():
            element = ""
            t = 0
            for key in sorted(self._AtomStep.atoms_info.keys()):
                if self._elements[key] != element:
                    element = self._elements[key]
                    t += 1
                self._AtomStep.atoms_info[key]["type"] = t
        with open(self._args.output, "w") as write_file:
            for positions in self._AtomStep.cartesian_position.transpose((1, 0, 2)):
                write_file.write(f"ITEM: TIMESTEP\n{self._step_counter}\n")
                write_file.write(f"ITEM: NUMBER OF ATOMS\n{len(self._elements)}\n")
                write_file.write("ITEM: BOX BOUNDS pp pp pp\n")
                write_file.write(f"0.00000000  {self._AtomStep.lattice[0][0]:.8f}\n0.00000000  {self._AtomStep.lattice[1][1]}\n0.00000000  {self._AtomStep.lattice[2][2]}\n")
                write_file.write("ITEM: ATOMS id mol type element q xu yu zu\n")
                for atom_idx, position in enumerate(positions):
                    write_file.write(f"{atom_idx+1} {self._AtomStep.atoms_info[atom_idx]['molecule']+1} {self._AtomStep.atoms_info[atom_idx]['type']} {self._elements[atom_idx]} {self._AtomStep.atoms_info[atom_idx]['q']:.4f} {position[0]:.8f} {position[1]:.8f} {position[2]:.8f}\n")
                self._step_counter += 1

class LAMMPS_Single_Point(Single_Point):
    def __init__(self, args=None, atoms=None, **kwargs):
        super().__init__(input_obj=args, atoms=atoms, **kwargs)
        
        if hasattr(self, "_args") and hasattr(self._args, "molecules"):
            self.read_molecules()
        if hasattr(args, "bridge"):
            self.bridge = args
    @property
    def bridge(self):
        pass
    @bridge.setter
    def bridge(self, software):
        self.args, self.title, self.elements, self.atom_step = software.args, software.title, software.elements, software.atom_step
    def read_elements(self):
        with open(self._args.elements) as read_file:
            for line in read_file:
                if "dump_modify" in line and "element" in line:
                    element_list = line.split()[3:]
                    break
        with open(self._args.input) as read_file:
            for line in read_file:
                if "atoms" in line:
                    self.__atom_sum = int(line.split()[0])
                elif "Atoms" in line:
                    line = next(read_file)
                    for _ in range(self.__atom_sum):
                        sp = next(read_file).split()
                        self._elements.append(element_list[int(sp[2])-1])
                    self._elements = tuple(self._elements)
                    break
        return self._elements
    def read_all(self, atoms="all", rearange=True):
        if atoms == "all":
            self.atoms.bridge = self
            self.atoms.put(atoms)
        else:
            self.atoms = atoms
        atoms = (atoms if atoms == "all" else self._AtomStep.atom_list)
        if atoms != "all":
            atoms_dict = {atom: atom in atoms for atom in range(len(self._elements))}
        self._AtomStep.cartesian_flag = True
        if rearange:
            with open(self._args.elements) as read_file:
                for line in read_file:
                    if "dump_modify" in line and "element" in line:
                        element_list = line.split()[3:]
                        break
            elements_dict, elements_mole_dict = defaultdict(list), defaultdict(list)
        with open(self._args.input) as read_file:
            for line in read_file:
                if "atoms" in line:
                    atom_sum = int(line.split()[0])
                elif "xlo" in line:
                    i = 0
                    lattice = zeros((3, 3))
                    sp = line.split()
                    CellLowerX, CellUpperX = float(sp[0]), float(sp[1])
                    lattice[i, i] = CellUpperX - CellLowerX
                    sp = next(read_file).split()
                    i += 1
                    CellLowerY, CellUpperY = float(sp[0]), float(sp[1])
                    lattice[i, i] = CellUpperY - CellLowerY
                    sp = next(read_file).split()
                    i += 1
                    CellLowerZ, CellUpperZ = float(sp[0]), float(sp[1])
                    lattice[i, i] = CellUpperZ - CellLowerZ
                    self._AtomStep.lattice = lattice
                elif "Atoms" in line:
                    line = next(read_file)
                    for atom in range(atom_sum):
                        sp = next(read_file).split()
                        if not (atoms == "all" or atoms_dict[int(sp[0])-1]):
                            continue
                        if rearange:
                            elements_dict[element_list[int(sp[2])-1]].append([int(sp[0])-1, float(sp[4])+CellLowerX, float(sp[5])+CellLowerY, float(sp[6])+CellLowerZ])
                            elements_mole_dict[element_list[int(sp[2])-1]].append(int(sp[1])-1)
                        else:
                            self._AtomStep.put([float(sp[4])+CellLowerX, float(sp[5])+CellLowerY, float(sp[6])+CellLowerZ], int(sp[0])-1)
                            self._AtomStep.molecules.molecule_dictionary[int(sp[1])-1].append(int(sp[0])-1)
                    if rearange:
                        num = 0
                        elements, self._elements = [], []
                        cart_pos = []
                        for key, values in elements_dict.items():
                            for idx, val in enumerate(values):
                                self._elements.append(key)
                                elements.append(self._AtomStep.atoms.elements[val[0]])
                                cart_pos.append(val[1:])
                                self._AtomStep.molecules.molecule_dictionary[elements_mole_dict[key][idx]].append(num)
                                num += 1
                        self._AtomStep.elements, self._AtomStep.cartesian_position = elements, array(cart_pos)
                        self._elements = tuple(self._elements)
                    else:
                        self._AtomStep.lock()
                    break
        with open("molecules.csv", "w") as write_file:
            for mol_idx in sorted(self._AtomStep.molecules.molecule_dictionary.keys()):
                write_file.write("".join([str(atom) if idx == 0 else f",{atom}" for idx, atom in enumerate(self._AtomStep.molecules.molecule_dictionary[mol_idx])])+'\n')
    def write_all(self):
        with open(self._args.molecules) as read_file:
            row = 0
            for line in read_file:
                if compile(r"^\d+(,\d)*$").match(line.replace('\n', '')) is not None:
                    self._AtomStep.molecules.molecule_dictionary[row] = [int(_) for _ in line.replace('\n', '').split(",")]
                    row += 1
        from ase.data import atomic_masses, atomic_numbers
        with open(self._args.output, "w") as write_file:
            element_idx_dict = {}
            element_idx = 0
            former_element = ""
            recorded_elements = []
            for element, atom_element in zip(self._elements, self._AtomStep.atoms.elements):
                if (atom_element != former_element) and atom_element not in recorded_elements:
                    element_idx += 1
                    former_element = atom_element
                    atom_masses = atomic_masses[atomic_numbers[element.split("_")[0]]]
                    element_idx_dict[atom_element] = [element_idx, atom_masses]
                    recorded_elements.append(atom_element)
            write_file.write("LAMMPS data file\n\n")
            write_file.write(f"{len(self._elements)} atoms\n")
            write_file.write(f"{len(element_idx_dict)} atom types\n\n")
            write_file.write(f"{0:.10f} {self._AtomStep.lattice[0][0]*self._scaling:.10f} xlo xhi\n")
            write_file.write(f"{0:.10f} {self._AtomStep.lattice[1][1]*self._scaling:.10f} ylo yhi\n")
            write_file.write(f"{0:.10f} {self._AtomStep.lattice[2][2]*self._scaling:.10f} zlo zhi\n\n")
            write_file.write("Masses\n\n")
            for idx, val in enumerate(element_idx_dict.values(), 1):
                write_file.write(f"{idx} {val[1]:.3f}\n")
            write_file.write("\nAtoms\n\n")
            for atom in self._AtomStep.atom_list:
                for mol_idx, mol_atoms in self._AtomStep.molecules.molecule_dictionary.items():
                    if atom in mol_atoms:
                        write_file.write(f"{atom+1:6d}{mol_idx+1:6d}{element_idx_dict[self._AtomStep.atoms.elements[atom]][0]:6d}  0.0000  {self._AtomStep.cartesian_position[atom][0]*self._scaling:18.10e}{self._AtomStep.cartesian_position[atom][1]*self._scaling:18.10e}{self._AtomStep.cartesian_position[atom][2]*self._scaling:18.10e}\n")
