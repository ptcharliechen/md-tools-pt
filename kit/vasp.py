from os import getcwd
from numpy import zeros, array
from kit.fundamental import *

class POSCAR(Single_Point):
    def __init__(self, args=None, atoms=None):
        Single_Point.__init__(self, args, atoms)
        if hasattr(self, "_args") and hasattr(self._args, "input") and ("POSCAR" in self._args.input or "CONTCAR" in self._args.input):
            self.read_elements()
        self._relax_flag = False
    @property
    def bridge(self):
        pass
    @bridge.setter
    def bridge(self, software):
        self.relax = software.relax if hasattr(software, 'relax') else None
        self.relax_border = software.relax_border if hasattr(software, 'relax_border') else 0
        self._relax_flag = False if self.relax is None or not self.relax_border else True
        self.args, self.title, self.atom_step = software.args, software.title, software.atom_step
    def read_elements(self):
        with open(self._args.input) as read_file:
            for row, line in enumerate(read_file, 1):
                if row == 6:
                    elements = line.split()
                elif row == 7:
                    numbers = [int(num) for num in line.split()]
                    self._AtomStep.atoms.elements = Convert.eleNum2elements(elements, numbers)
                    for _ in range(len(self._AtomStep.atoms.elements)):
                        self._AtomStep.atoms_info[_]["element"] = self._AtomStep.atoms.elements[_]
                    self._AtomStep.molecules.elements = self._AtomStep.atoms.elements
                    return self._AtomStep.atoms.elements
    def rearrange(self, elements, positions, num=None, typ='f'):
        from collections import defaultdict
        assert len(elements) == len(positions), "Lengths of elements and AtomStep should be equal."
        tmp = defaultdict(list)
        for element, position in zip(elements, positions):
            tmp[element].append(position)
        if 'f' in typ.lower():
            self._AtomStep.fractional_flag, self._AtomStep.cartesian_flag = True, False
        elif 'c' in typ.lower():
            self._AtomStep.cartesian_flag, self._AtomStep.fractional_flag = True, False
        self._AtomStep.atoms.elements = list([element for element in tmp.keys() for _ in range(len(tmp[element]))])
        if num is None:
            atom_sum = len(self._AtomStep.atoms.elements)
            num = range(atom_sum)
        for atom in num:
            self._AtomStep.atoms_info[atom]["element"] = self._AtomStep.atoms.elements[atom]
        positions = []
        atom = 0
        for element in tmp.keys():
            for position in tmp[element]:
                positions.append(position)
                self._AtomStep.atoms.put(atom)
                atom += 1
        self._AtomStep.lock()
        if self._AtomStep.fractional_flag:
            self._AtomStep.fractional_position = positions
        elif self._AtomStep.cartesian_flag:
            self._AtomStep.cartesian_position = positions
    def read_all(self, atoms="all"):
        if atoms == "all":
            self.atoms.bridge = self
            self.atoms.put(atoms)
        else:
            self.atoms = atoms
        atoms = (atoms if atoms == "all" else self._AtomStep.atom_list)
        with open(self._args.input) as read_file:
            row = 0
            for line in read_file:
                if "direct" in line.lower() or "cartesian" in line.lower() or "selective" in line.lower() or line.replace('\n', '') == "S":
                    if "selective" in line.lower() or line.replace('\n', '') == "S":
                        line = next(read_file)
                    if "direct" in line.lower():
                        self._AtomStep.fractional_flag, self._AtomStep.cartesian_flag = True, False
                        self._AtomStep.bridge = self
                        self._relaxBorder = 0
                    elif "cartesian" in line.lower():
                        self._AtomStep.cartesian_flag, self._AtomStep.fractional_flag = True, False
                        self._AtomStep.bridge = self
                        self._relaxBorder = [0, 0, 0]
                    for atom in range(atom_sum):
                        sp = next(read_file).split()
                        row += 1
                        if (atoms == "all" or atom in atoms):
                            self._AtomStep.put([float(_) for _ in sp[:3]], atom)
                        if self._AtomStep.fractional_flag:
                            if len(sp) > 3:
                                self._relax_flag = True
                                if sp[5] == "F" and float(sp[2]) > self._relax_border:
                                    self._relax_border = float(sp[2])
                        elif self._AtomStep.cartesian_flag:
                            if len(sp) > 3:
                                self._relax_flag = True
                                if sp[5] == "T" and float(sp[2]) > self._relax_border[2]:
                                    self._relax_border = [float(_) for _ in sp[:3]]
                    self._AtomStep.lock()
                    if self._AtomStep.cartesian_flag and self._relax_flag:
                        from numpy import linalg
                        self._relax_border = array(self._relax_border).dot(linalg.inv(self._lattice))[2]
                elif row < 9:
                    row += 1
                    self._title = line.replace('\n', '')
                    
                    line = next(read_file)
                    row += 1
                    scale = float(line.split()[0])
                    
                    arr = zeros((3, 3))
                    for _ in range(3):
                        line = next(read_file)
                        row += 1
                        arr[_] = [scale*self._scaling*float(_) for _ in line.split()]
                    self._AtomStep.lattice = arr
                    
                    line = next(read_file)
                    row += 1
                    elements = line.split()
                    line = next(read_file)
                    self._AtomStep.atoms.elements = Convert.eleNum2elements(elements, [int(_) for _ in line.split()])
                    row += 1
                    atom_sum = len(self._AtomStep.atoms.elements)
                    for _ in range(atom_sum):
                        self._AtomStep.atoms_info[_]["element"] = self._AtomStep.atoms.elements[_]
                    self._AtomStep.bridge = self
    def read_fast(self, atoms="all"):
        if atoms == "all":
            self.atoms.bridge = self
            self.atoms.put(atoms)
        else:
            self.atoms = atoms
        atoms = (atoms if atoms == "all" else self._AtomStep.atom_list)
        self._AtomStep.fast_flag = True
        with open(self._args.input) as read_file:
            for line in read_file:
                if "direct" in line.lower() or "cartesian" in line.lower():
                    if "direct" in line.lower():
                        self._AtomStep.fractional_flag, self._AtomStep.cartesian_flag = True, False
                        self._AtomStep.bridge = self
                        self._relax_border = 0
                    elif "cartesian" in line.lower():
                        self._AtomStep.cartesian_flag, self._AtomStep.fractional_flag = True, False
                        self._AtomStep.bridge = self
                        self._relax_border = [0, 0, 0]
                    for atom in range(atom_sum):
                        sp = next(read_file).split()
                        if (atoms == "all" or atom in atoms):
                            self._AtomStep.fast_append(atom, [_ for _ in sp])
                    self._AtomStep.fast_step()
                    if self._AtomStep.cartesian_flag and self._relax_flag:
                        from numpy import linalg
                        self._relax_border = array(self._relax_border).dot(linalg.inv(self._lattice))[2]
                else:
                    self._title = line.replace('\n', '')
                    
                    line = next(read_file)
                    scale = float(line.split()[0])
                    arr = zeros((3, 3))
                    
                    for _ in range(3):
                        line = next(read_file)
                        arr[_] = [scale*self._scaling*float(_) for _ in line.split()]
                    self._AtomStep.lattice = arr
                    
                    line = next(read_file)
                    elements = line.split()
                    line = next(read_file)
                    self._AtomStep.atoms.elements = Convert.eleNum2elements(elements, [int(_) for _ in line.split()])
                    atom_sum = len(self._AtomStep.atoms.elements)
                    for _ in range(atom_sum):
                        self._AtomStep.atoms_info[_]["element"] = self._AtomStep.atoms.elements[_]
                    self._AtomStep.bridge = self
                    
                    line = next(read_file)
                    if "selective" in line.lower():
                        line = next(read_file)
    def write_all(self, typ='f'):
        lattice = self._scaling*self._AtomStep.lattice
        elements = self._AtomStep.atoms.elements
        if typ.lower() == 'f':
            position = self._AtomStep.fractional_position
            frac_pos = position
        elif typ.lower() == 'c':
            position = self._AtomStep.cartesian_position
            frac_pos = self._AtomStep.fractional_position
        with open(self._args.output, "w") as write_file:
            write_file.write("{}\n".format(self._title.replace('\n', '')))
            write_file.write("1.0\n")
            for _ in range(3):
                write_file.write(f"    {lattice[_][0]:.10f}\t{lattice[_][1]:.10f}\t{lattice[_][2]:.10f}\n")
            eleNum, element_order = Convert.elements2eleNum(elements)
            for element in element_order:
                write_file.write(f"     {element}")
            write_file.write("\n")
            for element in element_order:
                write_file.write(f"     {eleNum[element]}")
            atom_sum = len(elements)
            write_file.write("\n")
            if self._relax_flag:
                write_file.write("Selective dynamics\n")
            if typ.lower() == 'f':
                write_file.write("Direct\n")
            elif typ.lower() == 'c':
                write_file.write("Cartesian\n")
            for idx in range(atom_sum):
                write_file.write(f"    {position[idx][0]:.8f}\t{position[idx][1]:.8f}\t{position[idx][2]:.8f}\t")
                if self._relax_flag and isinstance(self._relax_border, list):
                    for _ in range(len(self._relax_border[idx])):
                        write_file.write(f" {('T' if self._relax_border[idx][_] == 'T' else 'F')}")
                elif self._relax_flag and isinstance(self._relax_border, float):
                    write_file.write((" F F F" if frac_pos[idx][2] <= self._relax_border else " T T T"))
                write_file.write('\n')
    def write_fast(self, typ='f'):
        lattice = self._scaling*self._AtomStep.lattice
        elements = self._AtomStep.atoms.elements
        if len(self._AtomStep.fast_position) == 1:
            frac_pos = self._AtomStep.fast_position[0]
        else:
            frac_pos = self._AtomStep.fast_position
            
        with open(self._args.output, "w") as write_file:
            write_file.write("{}\n".format(self._title.replace('\n', '')))
            write_file.write("1.0\n")
            for _ in range(3):
                write_file.write(f"    {lattice[_][0]:10f}\t{lattice[_][1]:10f}\t{lattice[_][2]:10f}\n")
            eleNum, element_order = Convert.elements2eleNum(elements)
            for element in element_order:
                write_file.write(f"     {element}")
            write_file.write("\n")
            for element in element_order:
                write_file.write(f"     {eleNum[element]}")
            atom_sum = len(elements)
            write_file.write("\n")
            if len(frac_pos[0]) > 3:
                write_file.write("Selective dynamics\n")
            if typ == 'f':
                write_file.write("Direct\n")
            elif typ == 'c':
                write_file.write("Cartesian\n")
            for idx in range(atom_sum):
                write_file.write(f"    {frac_pos[idx][0]}\t{frac_pos[idx][1]}\t{frac_pos[idx][2]}\t")
                if len(frac_pos[idx]) > 3:
                    write_file.write(" ".join([f"{pos}" for pos in frac_pos[idx][3:]]))
                if self._relax_flag and isinstance(self._relax_border, list):
                    for _ in range(len(self._relax_border[idx])):
                        write_file.write(" {}".format(("T" if self._relax_border[idx][_] == "T" else "F")))
                elif self._relax_flag and isinstance(self._relax_border, float):
                    write_file.write((" F F F" if frac_pos[idx][2] <= self._relax_border else " T T T"))
                write_file.write('\n')

class XDATCAR(Position, Periodic, Trajectory):
    def __init__(self, args=None, steps=None, atoms=None):
        Position.__init__(self, args)
        Periodic.__init__(self, args)
        Trajectory.__init__(self, args, steps, atoms)
        if hasattr(self, "_args") and hasattr(self._args, "input") and "XDATCAR" in self._args.input:
            self.read_elements()
    @property
    def bridge(self):
        pass
    @bridge.setter
    def bridge(self, software):
        self.args, self.title, self._AtomStep = software.args, software.title, software._AtomStep
    def read_elements(self):
        from os import popen
        with popen(f"head {self._args.input}") as command:
            for row, line in enumerate(command, 1):
                if row == 6:
                    elements = line.split()
                elif row == 7:
                    numbers = [int(num) for num in line.split()]
                    self._AtomStep.atoms.elements = Convert.eleNum2elements(elements, numbers)
                    for _ in range(len(self._AtomStep.atoms.elements)):
                        self._AtomStep.atoms_info[_]["element"] = self._AtomStep.atoms.elements[_]
                    self._AtomStep.molecules.elements = self._AtomStep.atoms.elements
                    return self._AtomStep.atoms.elements
    def read_fast(self, steps="all", atoms="all"):
        if atoms == "all":
            self.atoms.bridge = self
            self.atoms.put("all")
        else:
            self.atoms = atoms
        atoms = (atoms if atoms == "all" else self._AtomStep.atom_list)
        if steps == "all":
            self.steps.bridge = self
            self.steps.put("all")
        else:
            self.steps = steps
            self.cutoff_step = max(steps.get())
        steps = (steps if steps == "all" else self._AtomStep.steps.get())
        self._AtomStep.fast_flag = True
        with open(self._args.input) as read_file:
            row = 1
            atom = -1
            for line in read_file:
                if "Direct configuration" in line:
                    self._stepCounter += 1
                    if not (steps == "all" or self._stepCounter in steps):
                        for _ in range(atom_sum):
                            line = next(read_file)
                    else:
                        for atom in range(atom_sum):
                            line = next(read_file)
                            if (atoms == "all" or atom in atoms):
                                self._AtomStep.fast_append(atom, [_ for _ in line.split()])
                        self._AtomStep.fast_step()
                else:
                    row = 0
                    while(1):
                        if row == 0:
                            if self._stepCounter == 0 and self._title == line:
                                self._AtomStep.lattice.NpT_flag = True
                            elif self._stepCounter == 0:
                                self._title = line
                            line = next(read_file)
                            row += 1
                        elif row == 1:
                            scale = float(line.split()[0])
                            if (steps == "all" or self._stepCounter in steps or self._stepCounter == 0):
                                arr = zeros((3, 3))
                                for i in range(3):
                                    line = next(read_file)
                                    arr[i] = [scale*self._scaling*float(_) for _ in line.split()]
                                    row += 1
                                self._AtomStep.lattice = arr
                            else:
                                for _ in range(3):
                                    line = next(read_file)
                                    row += 1
                            line = next(read_file)
                        elif row == 4 and self._stepCounter == 0:
                            elements = line.split()
                            line = next(read_file)
                            row += 1
                            self._AtomStep.atoms.elements = Convert.eleNum2elements(elements, [int(_) for _ in line.split()])
                            atom_sum = len(self._AtomStep.atoms.elements)
                            for _ in range(atom_sum):
                                self._AtomStep.atoms_info[_]["element"] = self._AtomStep.atoms.elements[_]
                            self._AtomStep.fractional_flag, self._AtomStep.cartesian_flag = True, False
                            self._AtomStep.bridge = self
                            break
    def read_all(self, steps="all", atoms="all"):
        if atoms == "all":
            self.atoms.bridge = self
            self.atoms.put(atoms)
        else:
            self.atoms = atoms
        atoms = (atoms if atoms == "all" else self._AtomStep.atom_list)
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
            StepsDict = {step: step in steps for step in range(1, self.steps.total_steps+1)}
        self._AtomStep.fractional_flag, self._AtomStep.cartesian_flag = True, False
        with open(self._args.input) as read_file:
            self._title = None
            for line in read_file:
                if "Direct configuration" in line:
                    self._stepCounter += 1
                    if self._stepCounter % 1000 == 0:
                        print(f"{self._stepCounter} / {self._cutoff_step if self._cutoff_step != -1 else self._AtomStep.steps.total_steps}")
                    if not (steps == "all" or StepsDict[self._stepCounter]):
                        for _ in range(atom_sum):
                            line = next(read_file)
                    else:
                        for atom in range(atom_sum):
                            line = next(read_file)
                            if (atoms == "all" or AtomsDict[atom]):
                                self._AtomStep.put([float(_) for _ in line.split()], atom)
                else:
                    row = 0
                    while(1):
                        if row == 0:
                            if self._stepCounter == 0 and self._title == line:
                                self._AtomStep.lattice.NpT_flag = True
                            elif self._stepCounter == 0:
                                self._title = line
                            line = next(read_file)
                            row += 1
                        elif row == 1:
                            scale = float(line.split()[0])
                            if (steps == "all" or self._stepCounter in steps or self._stepCounter == 0):
                                arr = zeros((3, 3))
                                for i in range(3):
                                    line = next(read_file)
                                    arr[i] = [scale*self._scaling*float(_) for _ in line.split()]
                                    row += 1
                                self._AtomStep.lattice = arr
                            else:
                                for _ in range(3):
                                    line = next(read_file)
                                    row += 1
                            line = next(read_file)
                        elif row == 4 and self._stepCounter == 0:
                            elements = line.split()
                            line = next(read_file)
                            row += 1
                            self._AtomStep.atoms.elements = Convert.eleNum2elements(elements, [int(_) for _ in line.split()])
                            atom_sum = len(self._AtomStep.atoms.elements)
                            for _ in range(atom_sum):
                                self._AtomStep.atoms_info[_]["element"] = self._AtomStep.atoms.elements[_]
                            if atoms != "all":
                                AtomsDict = {atom: atom in atoms for atom in range(atom_sum)}
                            self._AtomStep.fractional_flag, self._AtomStep.cartesian_flag = True, False
                            self._AtomStep.bridge = self
                            break
                if (self._stepCounter-1) > self._cutoff_step and self._cutoff_step != -1:
                    break
            self._AtomStep.fractional_flag = True
            self._AtomStep.lock()
    def write_all(self):
        self._lattice = self._scaling*self._AtomStep.lattice
        frac_pos = self._AtomStep.fractional_position.transpose((1, 0, 2))
        with open(self._args.output, "w") as write_file:
            flag = True
            for i in range(len(frac_pos)):
                if flag:
                    write_file.write("{}\n".format(self._title.replace('\n', '')))
                    write_file.write("1.0\n")
                    if self._lattice.ndim == 3:
                        for j in range(3):
                            write_file.write(f"    {self._lattice[i][j][0]:.10f}\t{self._lattice[i][j][1]:.10f}\t{self._lattice[i][j][2]:.10f}\n")
                    else:
                        for j in range(3):
                            write_file.write(f"    {self._lattice[j][0]:.10f}\t{self._lattice[j][1]:.10f}\t{self._lattice[j][2]:.10f}\n")
                    eleNum, element_order = Convert.elements2eleNum(self._AtomStep.atoms.elements)
                    for element in element_order:
                        write_file.write(f"{element:>8}")
                    write_file.write('\n')
                    for element in element_order:
                        write_file.write(f"{eleNum[element]:>8}")
                    write_file.write('\n')
                if self._lattice.ndim == 2:
                    flag = False
                write_file.write(f"Direct configuration=    {i+1}\n")
                for j in range(len(self._AtomStep.atoms.elements)):
                    write_file.write(f"    {frac_pos[i][j][0]:.8f}\t{frac_pos[i][j][1]:.8f}\t{frac_pos[i][j][2]:.8f}\n")
    def write_fast(self):
        frac_pos = self._AtomStep.fast_position
        with open(self._args.output, "w") as write_file:
            flag = True
            for i in range(len(frac_pos)):
                if flag:
                    write_file.write("{}\n".format(self._title.replace('\n', '')))
                    write_file.write("1.0\n")
                    if self._lattice.ndim == 3:
                        for j in range(3):
                            write_file.write(f"    {self._lattice[i][j][0]}\t{self._lattice[i][j][1]}\t{self._lattice[i][j][2]}\n")
                    else:
                        for j in range(3):
                            write_file.write(f"    {self._lattice[j][0]}\t{self._lattice[j][1]}\t{self._lattice[j][2]}\n")
                    eleNum, element_order = Convert.elements2eleNum(self._AtomStep.atoms.elements)
                    for element in element_order:
                        write_file.write(f"{element:>8}")
                    write_file.write('\n')
                    for element in element_order:
                        write_file.write(f"{eleNum[element]:>8}")
                    write_file.write('\n')
                if self._lattice.ndim == 2:
                    flag = False
                write_file.write(f"Direct configuration=    {i+1}\n")
                for j in range(len(self._AtomStep.atoms.elements)):
                    write_file.write(f"    {frac_pos[i][j][0]}\t{frac_pos[i][j][1]}\t{frac_pos[i][j][2]}\n")

class ACF(Fundamental, Charge):
    def __init__(self, input_obj, atoms=None):
        Fundamental.__init__(self, input_obj)
        Charge.__init__(self)
        if hasattr(input_obj, "bridge"):
            self.bridge = input_obj.bridge
        self._AtomStep = AtomStep_Single_Point(input_obj, atoms)
    @property
    def bridge(self):
        pass
    @property
    def diff(self):
        self.read_ref()
        ref = [self._ref[self._elements[atom]] for atom in self.atoms.get()]
        return self._charge - array(ref)
    @bridge.setter
    def bridge(self, software):
        self.args, self.elements, self.reference = software.args, software.elements, software.reference
    def read_elements(self):
        charge_edit = Charge_Edit(self._args)
        charge_edit.read_elements()
        self._elements = charge_edit.elements
    def read_ref(self):
        charge_edit = Charge_Edit(self._args)
        charge_edit.read_ref()
        self._ref = charge_edit.reference
    def read_all(self, file):
        if path.isfile(path.join(file, "ACF.dat")):
            file = path.join(file, "ACF.dat")
        with open(file) as read_file:
            row = 0
            for line in read_file:
                row += 1
                if row > 2 and row < len(self._AtomStep.atoms.elements)+3:
                    self._charge.append(-float(line.split()[4]))
        self._charge = array(self._charge)
    def read_atom(self, file):
        if path.isfile(path.join(file, "ACF.dat")):
            file = path.join(file, "ACF.dat")
        with open(file) as read_file:
            tmp = []
            row = 0
            for line in read_file:
                row += 1
                if row > 2 and row-3 in self.atoms.get():
                    tmp.append(-float(line.split()[4]))
        self._charge = tmp

class ACFs(Trajectory, Charge):
    def __init__(self, input_obj, steps=None, atoms=None):
        Trajectory.__init__(self, input_obj, steps=steps, atoms=atoms)
        Charge.__init__(self)
        if hasattr(input_obj, "bridge"):
            self.bridge = input_obj
            self.steps = Step(input_obj) if steps is None else steps
            self.atoms = Atom(input_obj) if atoms is None else atoms
    @property
    def bridge(self):
        pass
    @bridge.setter
    def bridge(self, software):
        self.args, self.atom_step, self.reference = software.args, software.atom_step, software.reference
    def read_elements(self):
        charge_edit = Charge_Edit(self._args)
        charge_edit.read_elements()
        self._AtomStep.elements = charge_edit.elements
    def read_ref(self):
        charge_edit = Charge_Edit(self._args)
        charge_edit.read_ref()
        self._ref = charge_edit.reference
    def read_all(self, diff_file="Charge_Diff"):
        from etc import read_csv
        if path.splitext(diff_file)[1] == "":
            diff_file += ".csv"
        self._diff, steps = read_csv(path.join(self._args.input, diff_file), num_flag=True)
        for step in steps:
            self._AtomStep.steps.put(step)
    def read_atom(self, atoms=None, diff_file="Charge_Diff"):
        diff = defaultdict(list)
        self.atoms = atoms
        atom_list = atoms.get()
        if path.splitext(diff_file)[1] == "":
            diff_file += ".csv"
        steps = []
        with open(path.join(self._args.input, diff_file)) as read_file:
            for row, line in enumerate(read_file, 1):
                if row > 1:
                    steps.append(int(line.split(',')[0]))
                    for atom, charge in enumerate(line.split(',')[1:]):
                        if atom in atom_list:
                            diff[atom].append(float(charge))
        self._AtomStep.steps.total_steps = max(steps)
        for step in steps:
            self._AtomStep.steps.put(step)
        for atom in atom_list:
            self._diff.append(diff[atom])
        self._diff = array(self._diff).T
    def write_all(self):
        atom_list = self.atoms.get()
        with open(path.join(getcwd(), f"{self._args.output}.csv"), "a+") as write_file:
            for atom in atom_list:
                write_file.write(f",{atom+self._args.atom}")
            write_file.write(",sum\n")
            write_file.write("ref")
            for atom in atom_list:
                write_file.write(f",{self._ref[self._AtomStep.elements[atom]]}")
            write_file.write(",0\n")
            for idx, step in enumerate(self._AtomStep.steps.get()):
                tmp = []
                write_file.write(f"{step}")
                for _ in range(len(atom_list)):
                    tmp.append(self._diff[idx][_])
                    write_file.write(f",{self._diff[idx][_]}")
                write_file.write(f",{sum(tmp):.4f}\n")
            write_file.write('\n')

class Charge_Edit(Trajectory, Charge):
    def __init__(self, input_obj, steps=None, atoms=None):
        Trajectory.__init__(self, input_obj, steps, atoms)
        Charge.__init__(self)
        self._ref = {}
        if hasattr(input_obj, "bridge"):
            self.bridge = input_obj
            if steps is not None:
                self.steps = steps
            if atoms is not None:
                self.atoms = atoms
    def bridge(self):
        pass
    def read_elements(self):
        if path.isfile(path.join(self._args.input)):
            dir_name = path.dirname(self._args.input)
        elif path.isdir(path.join(self._args.input)):
            dir_name = self._args.input
        if not path.isfile(path.join(dir_name, "Atoms.txt")):
            if path.isfile(self._args.poscar):
                pass
            elif path.isfile(path.join(self._args.poscar, "POSCAR")):
                self._args.poscar = path.join(self._args.poscar, "POSCAR")
            elif path.isfile(path.join(self._args.poscar, "CONTCAR")):
                self._args.poscar = path.join(self._args.poscar, "CONTCAR")
            elif path.isfile(path.join(self._args.poscar, "XDATCAR")):
                self._args.poscar = path.join(self._args.poscar, "XDATCAR")
            else:
                raise Exception("Provide information about elements or the numbers of them.")
            tmp = self._args
            tmp.input = self._args.poscar
            poscar = POSCAR(tmp)
            poscar.read_elements()
            with open(path.join(dir_name, "Atoms.txt"), "w") as write_file:
                for element in poscar.elements:
                    write_file.write(f"{element} ")
        read_file = open(path.join(dir_name, "Atoms.txt"))
        self._AtomStep.elements = read_file.readline().strip().split()
        read_file.close()
    def read_ref(self):
        if path.isfile(path.join(self._args.input)):
            dir_name = path.dirname(self._args.input)
        elif path.isdir(path.join(self._args.input)):
            dir_name = self._args.input
        if not path.isfile(path.join(dir_name, "Ref.txt")):
            if path.isfile(self._args.potcar):
                pass
            elif path.isfile(path.join(self._args.potcar, "POTCAR")):
                self._args.potcar = path.join(self._args.potcar, "POTCAR")
            else:
                raise Exception("Provide information about charge.")
            with open(path.join(dir_name, "Ref.txt"), "w") as write_file:
                with open(self._args.potcar) as read_file:
                    chargeFlag = False
                    for row, line in enumerate(read_file, 1):
                        if "End of Dataset" in line:
                            chargeFlag = True
                        elif chargeFlag and len(line.split()) == 3:
                            write_file.write(line.split()[1])
                        elif chargeFlag and len(line.split()) == 1:
                            write_file.write(f"\t{line.split()[0]}\n")
                            chargeFlag = False
                        elif row == 1:
                            write_file.write(line.split()[1])
                        elif row == 2:
                            write_file.write(f"\t{line.split()[0]}\n")
        with open(path.join(dir_name, "Ref.txt")) as read_file:
            for line in read_file:
                self._ref[line.split()[0]] = -float(line.split()[1])
    def build(self):
        self.read_elements()
        self.read_ref()
        if not path.isdir(self._args.input):
            raise Exception(f"{self._args.input} should be a directory.")
        from os import listdir
        charges = {}
        steps = []
        for file in listdir(self._args.input):
            if path.isfile(path.join(self._args.input, file)):
                continue
            if not file:
                print(f"Warning: {path.abspath(path.join(self._args.input, file))} does not use a positive integer as the directory name. Skip it.")
                continue
            if path.isfile(path.join(self._args.input, file, "ACF.dat")):
                acf = ACF(self._args)
                acf.elements = self._AtomStep.elements
                acf.read_all(path.join(self._args.input, file, "ACF.dat"))
                steps.append(int(file))
                charges[int(file)] = acf.charge
            else:
                print(f"Warning: {path.abspath(path.join(self._args.input, file))} does not have ACF.dat file. Skip it.")
        self._AtomStep.steps.total_steps = max(steps)
        for step in steps:
            self._AtomStep.steps.put(step)
        self._AtomStep.steps.sort()
        for step in self._AtomStep.steps.get():
            self._charge.append(charges[step])
        self.charge_diff()
    def charge_diff(self):
        ref = [self._ref[element] for element in self._AtomStep.elements]
        ref = [ref]*len(self._charge)
        self._diff = array(self._charge) - array(ref)
    def summary(self, charge="Charge", diff="Charge_Diff"):
        from kit.etc import write_csv
        if path.splitext(charge)[1] == "":
            charge += ".csv"
        if not path.isfile(path.join(self._args.input, charge)):
            columns = list(range(len(self._AtomStep.elements)))
            indice = ["Atom", "ref"] + list(self._AtomStep.steps.get())
            data = [[element for element in self._AtomStep.elements], [self._ref[element] for element in self._AtomStep.elements]]
            for _ in self._charge:
                data.append(_)
            write_csv(path.join(self._args.input, charge), data, indice, columns)
        if path.splitext(diff)[1] == '':
            diff += ".csv"
        if not path.isfile(path.join(self._args.input, diff)):
            write_csv(path.join(self._args.input, diff), self._diff, self._AtomStep.steps.get(), list(range(len(self._AtomStep.elements))))
