from numpy import zeros, array
import numpy as np
from os import path
from kit.fundamental import *
from kit.accelerate import *

class QE(Position, Periodic, Trajectory):
    def __init__(self, args):
        Fundamental.__init__(self)
        Position.__init__(self)
        Periodic.__init__(self)
        Trajectory.__init__(self, args)
        self.args_check(args)
    @property
    def bridge(self):
        pass
    @property
    def fractional_position(self):
        self._cartPos = array(self._cartPos)
        if not self.__Lattice.NpT_flag:
            return c2f(self._cartPos, self._lattice)
        else:
            return c2f_acc(self._cartPos, self._lattice)
    @bridge.setter
    def bridge(self, software):
        self.args, self.lattice, self.elements, self.cartesian_position = software.args, software.lattice, software.elements, software.cartesian_position
    @fractional_position.setter
    def fractional_position(self, fracPos):
        self._fracPos = array(fracPos)
    def read_elements(self):
        assert path.isfile(self._args.element), "Provide information about elements. They could be in 'QE.in'."
        with open(self._args.element) as read_file:
            elemFlag = False
            for line in read_file:
                if "ATOMIC_POSITIONS" in line:
                    elemFlag = True
                elif elemFlag and len(line.split()) < 4:
                    elemFlag = False
                elif elemFlag:
                    self._elements.append(line.split()[0])
            self.__atomSum = len(self._elements)
    def read_all(self):
        flag = False
        self.__Lattice = Lattice()
        arr = zeros((3, 3))
        self.read_elements()
        with open(self._args.input) as read_file:
            num = -4
            for line in read_file:
                if len(line.split()) == 0 and self._stepCounter > 0:
                    flag = False
                    self._cartPos.append(atomPos)
                    continue
                elif len(line.split()) == 0 and self._stepCounter == 0:
                    flag = False
                    continue
                elif not flag:
                    flag = True
                    num = -4
                    self._stepCounter += 1
                    if self.__Lattice.NpT_flag and self._stepCounter == 1:
                        arr = zeros((3, 3))
                    atomPos = zeros((self.__atomSum, 3))
                    continue
                num += 1
                if num > -1:
                    atomPos[num] = [self._scaling*float(_) for _ in line.split()]
                    continue
                elif num < 0:
                    arr[num+3] = [self._scaling*float(_) for _ in line.split()]
                if num == -1:
                    self.__Lattice.lattice = arr
                if self._stepCounter > self._cutoffStep and self._cutoffStep != -1:
                    break
            if len(line.split()) != 0:
                self._cartPos.append(atomPos)
        self._lattice = self.__Lattice.lattice
    def write_all(self):
        self._args.output = self._title
        if "mdtrj" not in self._args.output:
            self._args.output += ".mdtrj"
        with open(self._args.output, "w") as write_file:
            for idx, cartPos in enumerate(self._cartPos):
                write_file.write("{:.12f}\t{:.12f}\t{:.12f}\n".write(0, 0, 0))
                if self._lattice.ndim == 3:
                    for _ in self._lattice[idx]:
                        write_file.write(f"{_[0]:.12f}\t{_[1]:.12f}\t{_[2]:.12f}\n")
                else:
                    for _ in self._lattice:
                        write_file.write(f"{_[0]:.12f}\t{_[1]:.12f}\t{_[2]:.12f}\n")
                for position in cartPos:
                    write_file.write(f"{position[0]:.12f}\t{position[1]:.12f}\t{position[2]:.12f}\n")
                write_file.write('\n')

class QE_Single_Point(Single_Point, Fundamental):
    def __init__(self, args):
        Fundamental.__init__(self)
        Single_Point.__init__(self, args)
        self.args_check(args)
    @property
    def bridge(self):
        pass
    @bridge.setter
    def bridge(self, software):
        self.args, self.elements, self.lattice, self.relax, self.relax_border, self.fractional_position, self.cartesian_position = software.args, software.elements, software.lattice, software.relax, software.relax_border, software.fractional_position, software.cartesian_position
    def read_elements(self):
        with open(self._args.input) as read_file:
            eleFlag = False
            for line in read_file:
                if "ATOMIC_POSITIONS" in line.split():
                    eleFlag = True
                    continue
                elif len(line.split() == 0) and eleFlag:
                    break
                if eleFlag:
                    self._elements.append(line.split()[0])
    def read_all(self):
        self._relaxBorder = []
        with open(self._args.input) as read_file:
            latticeFlag, eleFlag, posFlag_cart, posFlag_frac = False, False, False, False
            for line in read_file:
                if "ATOMIC_POSITIONS" in line.split():
                    eleFlag = True
                    if "angstrom" in line.split() or "bohr" in line.split():
                        from numpy import linalg
                        posFlag_cart = True
                    elif "crystal" in line.split():
                        posFlag_frac = True
                    else:
                        raise Exception("Only access 'bohr', 'angstrom', or 'crystal' in 'ATOMIC_POSITIONS' row.")
                    continue
                elif "CELL_PARAMETERS" in line.split():
                    latticeFlag = True
                    latticeRow = 0
                    if "bohr" in line.split():
                        self._scaling = 0.529177249
                    continue
                elif len(line.split()) == 0 and eleFlag:
                    eleFlag = False
                    if posFlag_cart:
                        self._cartPos = array(self._cartPos)
                        posFlag_cart = False
                    elif posFlag_frac:
                        self._fracPos = array(self._fracPos)
                        posFlag_frac = False
                elif "prefix" in line.split():
                    self._title = line.replace("'", "").replace(",", "").replace("=", "").split()[1]
                if eleFlag:
                    sp = line.split()
                    self._elements.append(sp[0])
                    if posFlag_cart:
                        self._cartPos.append([float(_) for _ in sp[1:4]])
                    elif posFlag_frac:
                        self._fracPos.append([float(_) for _ in sp[1:4]])
                    if len(sp) > 4:
                        self._relaxFlag = True
                        relaxDirection = ""
                        for _ in range(4, 7):
                            relaxDirection += ("T" if sp[_] == "1" else "F")
                        self._relaxBorder.append(relaxDirection)
                    else:
                        self._relaxBorder.append("TTT")
                elif latticeFlag and latticeRow < 3:
                    self._lattice[latticeRow] = [float(_) for _ in line.split()]
                    latticeRow += 1
    def write_all(self, unit="crystal"):
        self._lattice = self._scaling*self._lattice
        if unit == "angstrom" or unit == "bohr":
            self._cartPos = self._scaling*self._cartPos
            position = self._cartPos
        elif unit == "crystal":
            position = self._fracPos
        with open(self._args.output, "w") as write_file:
            write_file.write(f"ATOMIC_POSITIONS {unit}\n")
            for idx, element in enumerate(self._elements):
                write_file.write(f"{element} {position[idx][0]:.8f} {position[idx][1]:.8f} {position[idx][2]:.8f}")
                if self._relaxFlag and isinstance(self._relaxBorder, list):
                    for _ in range(len(self._relaxBorder[idx])):
                        write_file.write(f" {'1' if self._relaxBorder[idx][_] == 'T' else '0'}")
                elif self._relaxFlag and isinstance(self._relaxBorder, float):
                    write_file.write((" 0 0 0" if position[idx][2] <= self._relaxBorder else " 1 1 1"))
                write_file.write('\n')
            write_file.write("\n\n")
            write_file.write(f"CELL_PARAMETERS {unit}\n")
            for _ in self._lattice:
                write_file.write(f"{_[0]:.10f}\t{_[1]:.10f}\t{_[2]:.10f}\n")

class Conquest(Position, Periodic, Trajectory):
    def __init__(self, args):
        Fundamental.__init__(self)
        Position.__init__(self)
        Periodic.__init__(self)
        Trajectory.__init__(self, args)
        self.args_check(args)
    @property
    def bridge(self):
        pass
    @property
    def fractional_position(self):
        self._cartPos = array(self._cartPos)
        if not self.__Lattice.NpT_flag:
            return c2f(self._cartPos, self._lattice)
        else:
            return c2f_acc(self._cartPos, self._lattice)
    @bridge.setter
    def bridge(self, software):
        self.args, self.lattice, self.elements, self.cartesian_position = software.args, software.lattice, software.elements, software.cartesian_position
    @fractional_position.setter
    def fractional_position(self, fracPos):
        self._fracPos = array(fracPos)
    def read_elements(self):
        eleFlag = False
        atomSum = -1
        from os import popen
        with popen(f"head {self._args.input}") as command:
            for line in command:
                if "PRIMCOORD" in line.upper():
                    eleFlag = True
                elif eleFlag and len(line.split()) == 2:
                    atomSum = int(line.split()[0])
                    break
        with popen(f"head -{atomSum+10} {self._args.input}") as command:
            for line in command:
                if "PRIMCOORD" in line.upper():
                    eleFlag = True
                elif eleFlag and atomCount == atomSum:
                    break
                elif eleFlag and len(line.split()) == 2:
                    atomCount = 0
                elif eleFlag and len(line.split()) == 4:
                    self._elements.append(line.split()[0])
                    atomCount += 1
    def read_all(self):
        atomPos = []
        self.__Lattice = Lattice()
        latticeFlag, posFlag = False, False
        num = -4
        with open(self._args.input) as read_file:
            for line in read_file:
                if "PRIMVEC" in line.upper():
                    latticeFlag = True
                    arr = zeros((3, 3))
                    self._stepCounter += 1
                elif "PRIMCOORD" in line.upper():
                    latticeFlag = False
                    self.__Lattice.lattice = arr
                    posFlag = True
                elif len(line.split()) == 0 and self._stepCounter > 0:
                    self.__atomSum = len(self._elements)
                    self._cartPos.append(atomPos)
                    atomPos = zeros((self.__atomSum, 3))
                    num = -4
                    posFlag = False
                elif num == -1 and posFlag:
                    num += 1
                elif self._stepCounter == 1 and posFlag:
                    self._elements.append(line.split()[0])
                    atomPos.append([float(_) for _ in line.split()[1:4]])
                elif posFlag:
                    atomPos[num] = [float(_) for _ in line.split()[1:]]
                    num += 1
                elif (latticeFlag and self.__Lattice.NpT_flag) or (latticeFlag and self._stepCounter == 1) or (latticeFlag and self._stepCounter == 2):
                    arr[num+4] = [self._scaling*float(_) for _ in line.split()]
                    num += 1
                elif latticeFlag:
                    num += 1
                if self._stepCounter > self._cutoffStep and self._cutoffStep != -1:
                    break
            if len(line.split()) != 0:
                self._cartPos.append(atomPos)
        self._lattice = self.__Lattice.lattice
    def write_all(self):
        with open(self._args.output, "w") as write_file:
            for idx, cartPos in enumerate(self._cartPos):
                write_file.write("CRYSTAL\n")
                write_file.write(f"PRIMVEC     \t{idx}\n")
                if self._lattice.ndim == 3:
                    for _ in self._lattice[idx]:
                        write_file.write(f"{_[0]:.12f}\t{_[1]:.12f}\t{_[2]:.12f}\n")
                else:
                    for _ in self._lattice:
                        write_file.write(f"{_[0]:.12f}\t{_[1]:.12f}\t{_[2]:.12f}\n")
                write_file.write(f"PRIMCOORD   \t{idx}\n")
                write_file.write(f"    {len(self._elements)}  1\n")
                for _, position in enumerate(cartPos):
                    write_file.write(f"  {self._elements[_]}\t{position[0]:.10f}\t{position[1]:.10f}\t{position[2]:.10f}\n")
                write_file.write('\n')

class Conquest_Single_Point(Single_Point, Fundamental):
    def __init__(self, args):
        Fundamental.__init__(self)
        Single_Point.__init__(self, args)
        self.args_check(args)
    @property
    def bridge(self):
        pass
    @bridge.setter
    def bridge(self, software):
        self.args, self.elements, self.lattice, self.relax, self.relax_border, self.fractional_position = software.preprocess, software.elements, software.lattice, software.relax, software.relax_border, software.fractional_position
    def read_elements(self):
        with open(self._args.element) as read_file:
            eleFlag = False
            elementOrder, elements = {}, {}
            for line in read_file:
                if "ChemicalSpeciesLabel" in line:
                    eleFlag = True
                    continue
                elif "endblock" in line and eleFlag:
                    break
                if eleFlag:
                    sp = line.split()
                    elementOrder[sp[0]] = sp[2]
                    elements[sp[2]] = 0
        with open(self._args.input) as read_file:
            row = 1
            for line in read_file:
                if row > 4:
                    sp = line.split()
                    elements[elementOrder[sp[3]]] += 1
        ele = [elementOrder[str(i)] for i in range(1, len(elementOrder.keys())+1)]
        num = [elements[elementOrder[str(i)]] for i in range(1, len(elementOrder.keys())+1)]
        self._elements = Convert.eleNum2elements(ele, num)
    def read_all(self):
        with open(self._args.element) as read_file:
            eleFlag = False
            elementOrder, elements = {}, {}
            for line in read_file:
                if "ChemicalSpeciesLabel" in line:
                    eleFlag = True
                    continue
                elif "endblock" in line and eleFlag:
                    break
                elif "IO.Title" in line.split():
                    if len(line.split()) == 1:
                        self._title = self._args.input.split('.')[0]
                    else:
                        self._title = line.split()[1]
                if eleFlag:
                    sp = line.split()
                    elementOrder[sp[0]] = sp[2]
                    elements[sp[2]] = 0
        with open(self._args.input) as read_file:
            row = 1
            self._relaxBorder = []
            for line in read_file:
                if row > 4 and row-5 < atomSum:
                    sp = line.split()
                    elements[elementOrder[sp[3]]] += 1
                    self._fracPos[row-5] = [float(_) for _ in line.split()[:3]]
                    if len(sp) > 3:
                        self._relaxFlag = True
                        relaxDirection = ""
                        for _ in range(4, 7):
                            relaxDirection += ("T" if sp[_].upper() == "T" else "F")
                        self._relaxBorder.append(relaxDirection)
                    else:
                        self._relaxBorder.append("TTT")
                elif 1 <= row < 4:
                    self._lattice[row-1] = [self._scaling*float(_) for _ in line.split()[:3]]
                elif row == 4:
                    atomSum = int(line.split()[0])
                    self._fracPos = zeros((int(line.split()[0]), 3))
                row += 1
        ele = [elementOrder[str(_)] for _ in range(1, len(elementOrder.keys())+1)]
        num = [elements[elementOrder[str(_)]] for _ in range(1, len(elementOrder.keys())+1)]
        self._elements = Convert.eleNum2elements(ele, num)
    def write_all(self):
        self._lattice = self._scaling*self._lattice
        with open(self._args.output, "w") as write_file:
            for _ in self._lattice:
                write_file.write(f"{_[0]:.10f}\t{_[1]:.10f}\t{_[2]:.10f}\n")
            write_file.write(f"{len(self._elements)}\n")
            element, element_num = self._elements[0], 1
            for idx, position in enumerate(self._fracPos):
                if idx > 0 and self._elements[idx-1] != self._elements[idx]:
                    element_num += 1
                write_file.write(f"{position[0]:.8f} {position[1]:.8f} {position[2]:.8f} {element_num} ")
                if self._relaxFlag and isinstance(self._relaxBorder, list):
                    for _ in range(len(self._relaxBorder[idx])):
                        write_file.write(f" {'T' if self._relaxBorder[idx][_] == 'T' else 'F'}")
                    write_file.write('\n')
                elif self._relaxFlag and isinstance(self._relaxBorder, float):
                    write_file.write((" F F F\n" if position[idx][2] <= self._relaxBorder else " T T T\n"))
                elif not self._relaxFlag:
                    write_file.write('T T T\n')

class Gromacs(Position, Periodic, Trajectory):
    def __init__(self, args):
        Position.__init__(self)
        Periodic.__init__(self)
        Trajectory.__init__(self, args)
        self.args_check(args)
    @property
    def bridge(self):
        pass
    @bridge.setter
    def bridge(self, software):
        self.args, self.elements, self.lattice, self.cartesian_position = software.args, software.elements, software.lattice, software.cartesian_position
    def read_elements(self):
        from os import popen
        with popen(f"grep -n 'TER' {self._args.input} | cut -d: -f1 | head -n 1") as command:
            row = command.read()
        with popen(f"head -{row} {self._args.input}") as command:
            for line in command:
                pass
    def read_all(self):
        from numpy import sin, cos, sqrt, deg2rad
        with open(self._args.input) as read_file:
            for line in read_file:
                if "CRYST" in line:
                    elements, cartPos, elementsIdx = [], [], []
                    elementsDict = {}
                    a, b, c, alpha, beta, gamma = [float(_) for _ in line.split()[1:7]]
                    alpha, beta, gamma = deg2rad(alpha), deg2rad(beta), deg2rad(gamma)
                    self._lattice = zeros((3, 3))
                    self._lattice[0][0] = a
                    self._lattice[1][0], self._lattice[1][1] = b*cos(deg2rad(alpha)), b*sin(deg2rad(alpha))
                    self._lattice[2][0], self._lattice[2][1],  = c*cos(deg2rad(beta)), c*(cos(gamma) - cos(alpha) * cos(beta))/sin(alpha)
                    self._lattice[2][2] = sqrt(c**2 - self._lattice[2][0]**2 - self._lattice[2][1]**2)
                elif "ATOM" in line:
                    elements.append(line.split()[-1])
                    cartPos.append([float(_) for _ in line.split()[-6:-3]])
                    if elements not in elementsIdx.keys():
                        elementsDict[elements] = 0
                    else:
                        elementsDict[elements] += 1
                    elementsIdx.append(elementsDict[elements])
                elif "TER" in line:
                    combine = list(zip(elements, elementsIdx, cartPos))
                    sortedCombine = sorted(combine, key=lambda x: x[0])
                    if self._stepCounter:
                        _, _, cartPos = zip(*sortedCombine)
                    else:
                        self._elements, _, cartPos = zip(*sortedCombine)
                    self._cartPos.append(cartPos)
                elif "MODEL" in line:
                    self._stepCounter += 1

class Gaussian(Single_Point):
    def __init__(self, args):
        Single_Point.__init__(self, args)
        self.__charge = 0
        self.__multiplicity = 1
        self.__chk = None
        self.__method = " cam-b3lyp/6-311g geom=connectivity"
        if args is not None:
            self.args_check(args)
    @property
    def bridge(self):
        pass
    @property
    def fractional_position(self):
        return c2f(self._AtomStep.cartesian_position, self._AtomStep.lattice)
    @bridge.setter
    def bridge(self, software):
        self.args, self.title, self.elements, self.atom_step = software.args, software.title, software.elements, software.atom_step
    @fractional_position.setter
    def fractional_position(self, fracPos):
        assert fracPos.ndim == 2, "Can not provide trajectory data."
        self._fracPos = array(fracPos)
    def read_elements(self):
        with open(self._args.input) as read_file:
            for line in read_file:
                if len(line.split()) == 2 and line.split()[0].isdigit() and line.split()[1].isdigit():
                    flag = True
                elif flag:
                    if len(line.split()) != 4:
                        flag = False
                    elif line.split()[0] == "Tv":
                        pass
                    else:
                        self._elements.append(line.split()[0])
    def read_all(self, atoms="all", rearange=True):
        from collections import defaultdict
        if atoms == "all":
            self.atoms.bridge = self
            self.atoms.put(atoms)
        else:
            self.atoms = atoms
        atoms = (atoms if atoms == "all" else self._AtomStep.atom_list)
        with open(self._args.input) as read_file:
            num = 0
            lattice, cartPos = [], []
            for line in read_file:
                if compile(r"-?\d+ \d+").match(line) is not None:
                    self.__charge = int(line.split()[0])
                    self.__multiplicity = int(line.split()[1])
                    line = next(read_file)
                    while(len(line.split()) != 4):
                        line = next(read_file)
                    if rearange:
                        elementsDict = defaultdict(list)
                    else:
                        elements = []
                    while(len(line.split()) == 4):
                        if line.split()[0] == "Tv":
                            lattice.append([float(_) for _ in line.split()[1:]])
                        elif num in atoms:
                            if rearange:
                                elementsDict[line.split()[0]].append([float(_) for _ in line.split()[1:]])
                            else:
                                elements.append(line.split()[0])
                                cartPos.append([float(_) for _ in line.split()[1:3]])
                        num += 1
                        line = next(read_file)
                    self._AtomStep.lattice = array(lattice)
                elif line[0] == '%':
                    self.__chk = line.split()[1:]
                elif line[0] == '#':
                    self.__method = line.split()[1:]
        num = 0
        if rearange:
            for key, values in elementsDict.items():
                for val in values:
                    self._elements.append(key)
                    cartPos.append(val)
                    num += 1
        self._AtomStep.elements, self._AtomStep.cartesian_position = self._elements, array(cartPos)
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
            for element, cartPos in zip(self._elements, self._AtomStep.cartesian_position):
                write_file.write(f"\t{element}\t{cartPos[0]:.8f}\t{cartPos[1]:.8f}\t{cartPos[2]:.8f}\n")
            if list(self._AtomStep.lattice) != []:
                for lattice in self._AtomStep.lattice:
                    write_file.write(f"\tTv\t{lattice[0]:10f}\t{lattice[1]:.10f}\t{lattice[2]:.10f}\n")

class LAMMPS(Position, Periodic, Trajectory):
    def __init__(self, args):
        Position.__init__(self)
        Periodic.__init__(self)
        Trajectory.__init__(self, args)
        self.args_check(args)
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
    def read_all(self, steps="all", atoms="all"):
        from os import popen
        if atoms == "all":
            self.atoms.bridge = self
            self.atoms.put(atoms)
        else:
            self.atoms = atoms
        atoms = (atoms if atoms == "all" else self._AtomStep.atom_list)
        if steps == "all":
            with popen(f"grep TIMESTEP {self._args.input} | wc") as command:
                self.steps.total_steps = int(command.read().split()[0])
            self.steps.put("all")
        else:
            self.steps = steps
            self.cutoff_step = max(steps.get())
        steps = (steps if steps == "all" else self._AtomStep.steps.get())
        self._AtomStep.cartesian_flag = True
        with open(self._args.input) as read_file:
            flag, stepFlag = True, False
            for line in read_file:
                if "TIMESTEP" in line:
                    line = next(read_file)
                    if self._stepCounter > self._cutoffStep and self._cutoffStep != -1:
                        break
                    self._stepCounter += 1
                    if not (steps == "all" or self._stepCounter in steps):
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
                    if not stepFlag:
                        atom_num = int(sp[0])
                        if atoms == "all":
                            self._AtomStep.atoms.put(tuple(range(atom_num)))
                        stepFlag = True
                elif "ATOMS" in line:
                    for atom in range(atom_num):
                        sp = next(read_file).split()
                        num = int(sp[0])-1
                        if (atoms == "all" or num in atoms):
                            self._AtomStep.put_step_atom_position(str(self._stepCounter), num, [np.float32(_) for _ in sp[-3:]], put_step_flag=False, put_atom_flag=False)
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
                write_file.write(f"ITEM: TIMESTEP\n{self._stepCounter}\n")
                write_file.write(f"ITEM: NUMBER OF ATOMS\n{len(self._elements)}\n")
                write_file.write("ITEM: BOX BOUNDS pp pp pp\n")
                write_file.write(f"0.00000000  {self._AtomStep.lattice[0][0]:.8f}\n0.00000000  {self._AtomStep.lattice[1][1]}\n0.00000000  {self._AtomStep.lattice[2][2]}\n")
                write_file.write("ITEM: ATOMS id mol type element q xu yu zu\n")
                for AtomIdx, position in enumerate(positions):
                    write_file.write(f"{AtomIdx+1} {self._AtomStep.atoms_info[AtomIdx]['molecule']+1} {self._AtomStep.atoms_info[AtomIdx]['type']} {self._elements[AtomIdx]} {self._AtomStep.atoms_info[AtomIdx]['q']:.4f} {position[0]:.8f} {position[1]:.8f} {position[2]:.8f}\n")
                self._stepCounter += 1
