from abc import ABC, abstractmethod
from collections.abc import Iterable
from os import path
from copy import deepcopy
from argparse import ArgumentParser, Namespace
from collections import defaultdict
from re import compile, split
from numpy import ndarray, array
import numpy as np

class Fundamental(ABC):
    def __init__(self, input_obj=None):
        self._title = "title"
        self._AtomStep = None
        if hasattr(input_obj, "bridge"):
            self.bridge = input_obj
        else:
            self.args = input_obj
    @property
    def args(self):
        return getattr(self, "_args", None)
    @property
    def steps(self):
        return self._AtomStep.steps
    @property
    def title(self):
        return self._title
    @abstractmethod
    def bridge(self):
        pass
    @property
    def elements(self):
        return self._AtomStep.atoms.elements
    @property
    def atoms_info(self):
        return self._AtomStep.atoms.atoms_info
    @property
    def atom_step(self):
        return self._AtomStep
    @property
    def lattice(self):
        return self._AtomStep.lattice
    @property
    def atoms(self):
        return self._AtomStep.atoms
    @property
    def molecules(self):
        return self._AtomStep.molecules
    @args.setter
    def args(self, args=None):
        self._args = deepcopy(Args.arg_check(args))
        if hasattr(self._args, "input"):
            self._title = self._args.input.split('.')[0] if hasattr(self._args, 'input') and self._args.input is not None else "position"
    @title.setter
    def title(self, title):
        self._title = str(title)
    @elements.setter
    def elements(self, elements):
        self._AtomStep.molecules.elements = self._AtomStep.atoms.elements = elements
    @atoms_info.setter
    def atoms_info(self, atoms_info):
        self._AtomStep.atoms.atoms_info = atoms_info
    @atom_step.setter
    def atom_step(self, atom_step):
        if atom_step.__class__.__bases__[0] == AtomStep:
            self._AtomStep = atom_step
        else:
            raise ValueError("Only 'AtomStep_Trj' or 'AtomStep_Single_Point' can be imported.")
    @lattice.setter
    def lattice(self, lattice):
        self._AtomStep.lattice = lattice
    @atoms.setter
    def atoms(self, atom_list):
        self._AtomStep.atoms = atom_list
    @molecules.setter
    def molecules(self, molecules):
        self._AtomStep.molecules = molecules
    def read_molecules(self, filepath=None):
        if filepath is None and not hasattr(self._args, "molecules"):
            raise ValueError("'molecules' argument does not exist.")
        elif filepath is None and self._args.molecules is None:
            while(1):
                filepath = input("Input the molecules file path: ")
                if path.isfile(filepath):
                    with open(filepath) as read_file:
                        for row in read_file:
                            if compile(r"^mol [a-zA-Z0-9]+$").match(row) is None and compile(r"^\d+(\,\d+)*$").match(row) is None:
                                print("Warning: Format of the molecules file is wrong.")
                                break
                        else:
                            self._args.molecules = filepath
                            break
                else:
                    print("Warning: File not found.")
        elif filepath is None:
            filepath = self._args.molecules
        
        self._AtomStep.molecules.read_molecules(filepath)
        
        self.atoms.molecule_kind = self._AtomStep.molecules.molecule_kind
        self.atoms.molecule_dictionary = self._AtomStep.molecules.molecule_dictionary
        for atom_idx, molecule in self._AtomStep.molecules.molecule_dictionary.items():
            for atom in molecule:
                self._AtomStep.atoms.atoms_info[atom]["molecule"] = atom_idx

class Position(Fundamental):
    def __init__(self, input_obj=None):
        Fundamental.__init__(self, input_obj)
        self._scaling = 1
    @property
    def unit_conversion(self):
        return self._scaling
    @property
    def cartesian_position(self):
        return self._AtomStep.cartesian_position
    @property
    def fast_position(self):
        return self._AtomStep.fast_position
    @unit_conversion.setter
    def unit_conversion(self, scaling):
        if compile(r'^\d+(\.\d+)?$').match(str(scaling)) is not None:
            self._scaling = float(scaling)
        else:
            raise ValueError("The scaling number should be a positive float.")
    @cartesian_position.setter
    def cartesian_position(self, cart_pos):
        self._AtomStep.cartesian_position = array(cart_pos)
    @fast_position.setter
    def fast_position(self, fast_pos):
        self._AtomStep.fast_position = array(fast_pos)

class Periodic(Fundamental):
    def __init__(self, input_obj=None):
        Fundamental.__init__(self, input_obj)
        from kit.accelerate import image_shift
        self.image_shift = image_shift
    @property
    def lattice(self):
        return self._AtomStep.lattice
    @property
    def fractional_position(self):
        return self._AtomStep.fractional_position
    @lattice.setter
    def lattice(self, lattice):
        self._AtomStep.lattice = lattice
    @fractional_position.setter
    def fractional_position(self, frac_pos):
        self._AtomStep.fractional_position = array(frac_pos)

class Trajectory(Fundamental):
    def __init__(self, input_obj=None, steps=None, atoms=None):
        self._stepCounter = 0
        self._cutoff_step = -1
        Fundamental.__init__(self, input_obj)
        self._AtomStep = AtomStep_Trj(input_obj, steps, atoms)
        self._title = "trajectory"
    @property
    def atom_step(self):
        return self._AtomStep
    @property
    def steps(self):
        return self._AtomStep.steps
    @property
    def steps_info(self):
        return self._AtomStep.steps_info
    @property
    def cutoff_step(self):
        pass
    @property
    def final_step(self):
        return self._stepCounter + self._args.step - 1
    @atom_step.setter
    def atom_step(self, AtomStep):
        if isinstance(AtomStep, AtomStep_Trj):
            self._AtomStep = AtomStep
        elif isinstance(AtomStep, AtomStep_Single_Point):
            self._AtomStep = AtomStep_Trj()
            self._AtomStep.bridge = AtomStep
            self._AtomStep.lattice = AtomStep.lattice
            self._AtomStep.steps.put("1")
            self._AtomStep.atoms = AtomStep.atoms
            if AtomStep.cartesian_flag:
                self._AtomStep.cartesian_position = AtomStep.cartesian_position
            if AtomStep.fractional_flag:
                self._AtomStep.fractional_position = AtomStep.fractional_position
            if AtomStep.fast_flag:
                self._AtomStep.fast_position = AtomStep.fast_position
        else:
            raise ValueError("Only 'AtomStep_Single_Point' or 'AtomStep_Trj' can be imported.")
    @steps.setter
    def steps(self, steps):
        self._AtomStep.steps = steps
    @cutoff_step.setter
    def cutoff_step(self, cutoff_step):
        if str(cutoff_step).isdigit() is not None:
            self._cutoff_step = int(cutoff_step) - self._args.step + 1
        else:
            raise ValueError("The cutoff step number should be a positive integer.")
    @steps_info.setter
    def steps_info(self, steps_info):
        self._AtomStep.steps_info = steps_info

class Single_Point(Periodic, Position):
    def __init__(self, input_obj=None, atoms=None):
        self._relax_flag = False
        self._relax_border = 0
        Periodic.__init__(self, input_obj)
        Position.__init__(self, input_obj)
        self._AtomStep = AtomStep_Single_Point(input_obj, atoms)
    @property
    def atom_step(self):
        return self._AtomStep
    @property
    def relax(self):
        return self._relax_flag
    @property
    def relax_border(self):
        return self._relax_border
    @atom_step.setter
    def atom_step(self, AtomStep):
        if isinstance(AtomStep, AtomStep_Single_Point):
            self._AtomStep = AtomStep
        elif isinstance(AtomStep, AtomStep_Trj):
            self._AtomStep = AtomStep_Single_Point()
            self._AtomStep.bridge = AtomStep
            self._AtomStep.lattice = AtomStep.lattice
            self._AtomStep.atoms = AtomStep.atoms
            if AtomStep.cartesian_flag:
                self._AtomStep.cartesian_position = AtomStep.cartesian_position[:, 0, :]
            if AtomStep.fractional_flag:
                self._AtomStep.fractional_position = AtomStep.fractional_position[:, 0, :]
            if AtomStep.fast_flag:
                self._AtomStep.fast_position = AtomStep.fast_position[:, 0, :]
        else:
            raise ValueError("Only 'AtomStep_Single_Point' or 'AtomStep_Trj' can be imported.")
    @relax.setter
    def relax(self, relax_flag):
        self._relax_flag = relax_flag
    @relax_border.setter
    def relax_border(self, relax_border):
        if compile(r"^-?\d+(\.\d+)?$").match(str(relax_border)) is not None:
            relax_border = float(relax_border)
            if -1 <= relax_border <= 1:
                self._relax_border = relax_border
            else:
                raise ValueError("The relax border should be between -1 and 1")
        else:
            raise ValueError("The relax border should give a number between -1 and 1")

class Charge:
    def __init__(self):
        self._charge, self._diff, self._ref = [], [], {}
    @property
    def reference(self):
        return self._ref
    @property
    def charge(self):
        return array(self._charge)
    @property
    def diff(self):
        return array(self._diff)
    @reference.setter
    def reference(self, ref):
        self._ref = ref

class Step:
    def __init__(self, input_obj=None, put_flag=False):
        self.__steps = []
        self._steps_info = defaultdict(dict)
        self._total_steps = -1
        self._index_dict = {}
        if hasattr(input_obj, "bridge"):
            self.bridge = input_obj
            if put_flag:
                if isinstance(input_obj, Step):
                    self.put(input_obj)
                else:
                    self.put(input_obj.steps)
        elif isinstance(input_obj, Args) or isinstance(input_obj, Namespace) or isinstance(input_obj, dict):
            self._args = Args.arg_check(input_obj)
            self._start = self._args.step
        else:
            self._start = 1
    @property
    def args(self):
        return getattr(self, "_args", None)
    @property
    def bridge(self):
        pass
    @property
    def start(self):
        return self._start
    @property
    def total_steps(self):
        return self._total_steps
    @property
    def steps_info(self):
        return self._steps_info
    @property
    def index_list(self):
        if self._index_dict == {}:
            self._index_dict = {step: index for index, step in enumerate(self.__steps)}
        return self._index_dict
    @start.setter
    def start(self, start):
        if compile(r'^\-?\d+$').match(str(start)) is not None:
            self._start = int(start)
        else:
            raise ValueError("The start number should be an integer.")
    @total_steps.setter
    def total_steps(self, total_steps):
        if str(total_steps).isdigit() is not None:
            self._total_steps = int(total_steps)
        else:
            raise ValueError("The total steps should be an integer.")
    @args.setter
    def args(self, args):
        self._args = Args.arg_check(args)
        if hasattr(self._args, "step"):
            self.start = self._args.step
        if self._total_steps == -1 and hasattr(self._args, "input"):
            self._set_total_steps()
    def _set_total_steps(self):
        from os import popen
        with popen(f"grep Direct {self._args.input} 2> /dev/null | wc") as command:
            self._total_steps = int(command.read().split()[0])
    @bridge.setter
    def bridge(self, software):
        if hasattr(software, "steps"):
            self.total_steps = software.steps.total_steps
        elif isinstance(software, Step):
            self.total_steps = software.total_steps
        self.args = software.args
    @steps_info.setter
    def steps_info(self, steps_info):
        if isinstance(steps_info, dict):
            self._steps_info = steps_info
        else:
            raise ValueError("The steps info should be 'dict'.")
    @index_list.setter
    def index_list(self, index_list):
        if isinstance(index_list, dict):
            self._index_dict = index_list
        else:
            raise ValueError("'index_list' should be 'dict'.")
    def put(self, input_steps):
        if compile(r"^\-?\d+$").match(str(input_steps)) is not None:
            input_steps = int(input_steps) - self._start + 1
            if input_steps > self._total_steps:
                return -3
            self._steps(input_steps)
        elif compile(r"^\-?\d+:\-?\d+(:\-?\d+)?$").match(str(input_steps)) is not None:
            sp = [int(step)-self._start+1 for step in input_steps.split(':')]
            if (len(sp) == 3 and int(sp[0]) > int(sp[2])) or (len(sp) == 2 and int(sp[0]) > int(sp[1])):
                return -2
            elif min(sp) < 0 and min(sp)+self._total_steps+1 < 0:
                return -4
            elif max(sp) > self._total_steps:
                return -3
            self._interval(input_steps)
        elif isinstance(input_steps, slice):
            self.__steps.append(input_steps)
        elif isinstance(input_steps, Iterable) and not isinstance(input_steps, str):
            if isinstance(input_steps[0], str):
                input_steps = [int(_) for _ in input_steps]
            elif not isinstance(input_steps[0], Iterable):
                for step in input_steps:
                    if step < 0:
                        self.__steps.append(step+1+self._total_steps)
                    else:
                        self.__steps.append(step)
            else:
                raise ValueError("Dimension should be one.")
            return 0
        elif compile(r"^[a-zA-z]+$").match(str(input_steps)) is not None:
            if input_steps.lower() == "end":
                return 1
            elif input_steps.lower() == "all":
                self.__steps = list(range(1, self._total_steps+1))
                return 2
            else:
                return -1
        else:
            return -1
        return 0
    def remove(self, input_steps):
        steps = Step(self)
        steps.put(input_steps)
        for step in steps.get():
            if step in self.__steps:
                self.__steps.remove(step)
    def get(self, start_flag=False, slice_flatten=False):
        if slice_flatten:
            steps = []
            for step in self.__steps:
                if isinstance(step, slice):
                    steps.extend(list(range(step.start, step.stop, step.step if step.step is not None else 1)))
                elif str(step).isdigit():
                    steps.append(step)
        else:
            steps = self.__steps
        if start_flag and self._start != 1:
            return tuple([step+self._start-1 for step in steps])
        else:
            return tuple(steps)
    def delete(self):
        self.__steps = []
    def sort(self):
        if self._start == 1:
            self.__steps = sorted(self.__steps)
        else:
            self.__steps = sorted([step-self._start+1 for step in self.__steps])
    @property
    def err_list(self):
        return defaultdict(int, {-4: "Warning: The negative numbers are too small.", -3: "Warning: Above the total steps.", -2: "Warning: The start number should be less than the end number.", -1: "Warning: Input error."})
    def _steps(self, steps):
        if isinstance(steps, Iterable):
            for step in steps:
                if step < 0:
                    self.__steps.append(step+self._total_steps-self._start+2)
                else:
                    self.__steps.append(step-self._start+1)
        elif isinstance(steps, int):
            if steps > 0:
                self.__steps.append(steps-self._start+1)
            else:
                self.__steps.append(steps+self._total_steps-self._start+2)
    def _interval(self, steps):
        step_range = steps.split(':')
        if int(step_range[0]) < 0:
            if len(step_range) == 2:
                step_range[0], step_range[1] = int(step_range[0])-self._start+self._total_steps+2, int(step_range[1])-self._start+self._total_steps+3
            else:
                step_range[0], step_range[1], step_range[2] = int(step_range[0])-self._start+self._total_steps+2, int(step_range[1]), int(step_range[2])-self._start+self._total_steps+3
        else:
            if len(step_range) == 2:
                step_range[0], step_range[1] = int(step_range[0])-self._start+1, int(step_range[1])-self._start+2
            else:
                step_range[0], step_range[1], step_range[2] = int(step_range[0])-self._start+1, int(step_range[1]), int(step_range[2])-self._start+2
        self.__steps.extend(list(range(step_range[0], step_range[2], step_range[1])) if len(step_range) == 3 else list(range(step_range[0], step_range[1])))

class Step_LAMMPS(Step):
    def __init__(self, input_obj=None, put_flag=False):
        Step.__init__(self, input_obj, put_flag)
    def _set_total_steps(self):
        from os import popen
        with popen(f"grep TIMESTEP {self._args.input} 2> /dev/null | wc") as command:
            self._total_steps = int(command.read().split()[0])

class Atom(Fundamental):
    def __init__(self, input_obj=None, put_flag=False):
        Fundamental.__init__(self, input_obj)
        del self._AtomStep
        self._atom_list = []
        self._index_dict = {}
        self._elements = []
        self._atoms_info = defaultdict(dict)
        self._molecules = Molecule(input_obj)
        if hasattr(input_obj, "bridge"):
            self.bridge = input_obj
            if put_flag:
                if isinstance(input_obj, Atom):
                    self.put(input_obj)
                else:
                    self.put(input_obj.atoms)
        elif isinstance(input_obj, Args) or isinstance(input_obj, Namespace) or isinstance(input_obj, dict):
            self._start = self._args.atom
        else:
            self._start = 0
        if hasattr(input_obj, "atoms_info"):
            self.atoms_info = input_obj.atoms_info
    @property
    def start(self):
        return self._start
    @property
    def bridge(self):
        pass
    @property
    def molecules(self):
        pass
    @property
    def elements(self):
        return self._elements
    @property
    def atoms_info(self):
        return self._atoms_info
    @property
    def index_list(self):
        if self._index_dict == {}:
            self._index_dict = {atom: index for index, atom in enumerate(self._atom_list)}
        return self._index_dict
    @property
    def molecule_kind(self):
        pass
    @property
    def molecule_dictionary(self):
        pass
    def atoms(self):
        raise ValueError("atoms can not be defined.")
    @bridge.setter
    def bridge(self, software):
        self.args, self.start = software.args, software.args.atom
        if hasattr(software, "atoms_info"):
            self.atoms_info = software.atoms_info
        if hasattr(software, "molecules"):
            self.molecules = software.molecules
        if software.elements == [] and hasattr(software, 'read_elements'):
            self.read_elements(software)
        elif software.elements != []:
            self.elements = software.elements
    @start.setter
    def start(self, start):
        if compile(r'^\-?\d+$').match(str(start)) is not None:
            self._start = int(start)
        else:
            raise ValueError("The start number should be an integer.")
    @elements.setter
    def elements(self, elements):
        if isinstance(elements, list):
            self._elements = elements
        elif isinstance(elements, dict):
            self._elements = Convert.eleNum2elements(elements.keys(), elements.values())
        else:
            raise ValueError("Elements should be 'list' or 'dict'.")
        self._atom_num = len(self._elements)
    @atoms_info.setter
    def atoms_info(self, atoms_info):
        if isinstance(atoms_info, dict):
            self._atoms_info = atoms_info
        else:
            raise ValueError("'atoms_info' should be 'dict'.")
    @index_list.setter
    def index_list(self, index_list):
        if isinstance(index_list, dict):
            self._index_dict = index_list
        else:
            raise ValueError("'index_list' should be 'dict'.")
    @molecules.setter
    def molecules(self, molecules):
        if isinstance(molecules, Molecule):
            self._molecules = molecules
        else:
            raise ValueError("'molecule' should be 'Molecule' object.")
    @molecule_kind.setter
    def molecule_kind(self, molecule_kind):
        self._molecules.molecule_kind = molecule_kind
    @molecule_dictionary.setter
    def molecule_dictionary(self, molecule_dictionary):
        self._molecules.molecule_dictionary = molecule_dictionary
    def put(self, input_atom_list):
        if compile(r"^\-?\d+$").match(str(input_atom_list)) is not None:
            if int(input_atom_list) > self._atom_num:
                return -3
            self._atoms(input_atom_list)
        elif isinstance(input_atom_list, str) and input_atom_list in self._elements:
            self._atom_kind(input_atom_list)
        elif compile(r"^\-?\d+:\-?\d+(:\-?\d+)?$").match(str(input_atom_list)) is not None:
            sp = [int(atom)-self._start for atom in input_atom_list.split(':')]
            if (len(sp) == 3 and int(sp[0]) > int(sp[2])) or (len(sp) == 2 and int(sp[0]) > int(sp[1])):
                return -2
            if min(sp) < 0 and min(sp)+self._atom_num < 0:
                return -4
            elif max(sp) > self._atom_num:
                return -3
            self._interval(input_atom_list)
        elif compile(r"^[a-zA-z0-9]+$").match(str(input_atom_list)) is not None:
            if input_atom_list.lower() == "end":
                return 1
            elif input_atom_list.lower() == "all":
                self._atom_list = list(range(self._atom_num))
                return 2
            elif "mol" in str(input_atom_list).lower():
                if self._molecules.molecule_dictionary == {}:
                    if not hasattr(self._args, "molecules"):
                        return -5
                    elif self._args.molecules is not None and path.isfile(self._args.molecules):
                        self._molecules.read_molecules(self._args.molecules)
                        self.read_atoms_info(self._args.molecules)
                    else:
                        return -6
                return self._mol(input_atom_list)
            else:
                return -1
        elif isinstance(input_atom_list, Iterable):
            if compile(r"^$").match(str(input_atom_list)) is not None:
                return -1
            elif compile(r"^\-?\d+$").match(str(input_atom_list[0])) is not None:
                self._atoms(input_atom_list)
            elif isinstance(input_atom_list, str):
                return -1
            else:
                raise ValueError("Dimension should be one.")
        elif isinstance(input_atom_list, slice):
            self._interval(f"{input_atom_list.start}:{input_atom_list.step if input_atom_list.step is not None else 1}:{input_atom_list.stop}")
        elif isinstance(input_atom_list, Atom):
            input_atom_list = input_atom_list.get()
            for atom in input_atom_list:
                self._atom_list.append(atom)
        else:
            return -1
        return 0
    def remove(self, input_atom_list):
        atoms = Atom(self)
        atoms.put(input_atom_list)
        for atom in atoms.get():
            if atom in self._atom_list:
                self._atom_list.remove(atom)
    def get(self, start_flag=False):
        atoms = self._atom_list
        if start_flag and self._start:
            return tuple([atom+self._start for atom in atoms])
        else:
            return tuple(atoms)
    def delete(self):
        self._atom_list = []
    def sort(self):
        self._atom_list = sorted(self._atom_list)
    @property
    def err_list(self):
        return defaultdict(int, {-7: "Warning: Can not match the molecule type.", \
                -6: "Warning: The molecule file does not exist.", \
                -5: "Warning: 'molecules' argument is not defined.", \
                -4: "Warning: The negative numbers are too small.", -3: "Warning: Above the total steps.", \
                -2: "Warning: The start number should be less than the end number.", -1: "Warning: Input error."})
    def _atoms(self, atoms):
        if isinstance(atoms, int) or isinstance(atoms, str):
            atoms = int(atoms)
            if atoms >= 0:
                self._in_list_check(atoms-self._start)
            else:
                self._in_list_check(atoms+self._atom_num-self._start)
        elif isinstance(atoms, Iterable):
            for atom in atoms:
                if int(atom) >= 0:
                    self._in_list_check(int(atom)-self._start)
                else:
                    self._in_list_check(int(atom)+self._atom_num-self._start)
    def _mol(self, inp):
        inp = inp.replace("mol", "").replace("_", "", 1)
        if inp in [_.split('_')[0] for _ in self._molecules.molecule_kind]:
            for i, kind in enumerate(self._molecules.molecule_kind):
                if kind.split('_')[0] == inp:
                    for atom in self._molecules.molecule_dictionary[i]:
                        self._in_list_check(atom)
            return 0
        elif inp in self._molecules.molecule_kind:
            for atom in self._molecules.molecule_dictionary[self._molecules.molecule_kind.index(inp)]:
                self._in_list_check(atom)
            return 0
        elif len(inp.split('_')) == 2 and inp.split('_')[1] in self._elements:
            for mol, kind in enumerate(self._molecules.molecule_kind):
                if inp.split('_')[0] == kind.split('_')[0]:
                    for atom in self._molecules.molecule_dictionary[mol]:
                        if inp.split('_')[1] == self._elements[atom]:
                            self._in_list_check(atom)
            return 0
        else:
            return -7
    def _interval(self, inp):
        atom_range = inp.split(':')
        if int(atom_range[0]) < 0:
            if len(atom_range) == 2:
                atom_range[0], atom_range[1] = int(atom_range[0])+self._atom_num-self._start, int(atom_range[1])+self._atom_num-self._start+1
            elif len(atom_range) == 3:
                atom_range[0], atom_range[1], atom_range[2] = int(atom_range[0])+self._atom_num-self._start, int(atom_range[1]), int(atom_range[2])+self._atom_num-self._start+1
            else:
                print("Warning: Two or three numbers separated by ':'.")
        else:
            if len(atom_range) == 2:
                atom_range[0], atom_range[1] = int(atom_range[0])-self._start, int(atom_range[1])-self._start+1
            elif len(atom_range) == 3:
                atom_range[0], atom_range[1], atom_range[2] = int(atom_range[0])-self._start, int(atom_range[1]), int(atom_range[2])-self._start+1
            else:
                print("Warning: Two or three numbers separated by ':'.")
        self._atom_list.extend(list(range(atom_range[0], atom_range[1])) if len(atom_range) == 2 else list(range(atom_range[0], atom_range[2], atom_range[1])))
    def _atom_kind(self, inp):
        element = inp.split('_')[0]
        tmp = []
        for idx, atom in enumerate(self._elements):
            if element == atom:
                tmp.append(idx)
        if inp in self._elements:
            for atom in tmp:
                self._in_list_check(atom)
        elif inp.split('_')[1].lower() == "except":
            except_num = []
            for split in inp.split('_')[2:]:
                if split.isdigit():
                    except_num += [int(split)-self._start]
                elif compile(r"^\d+(:\d+){1,2}$").match(split):
                    sp = split.split(':')
                    except_num += list(range(int(sp[0])-self._start, int(sp[2])-self._start+1, int(sp[1]))) if len(sp) == 3 else list(range(int(sp[0])-self._start, int(sp[1])-self._start+1))
                elif compile(r"^-\d+(:-?\d+)?(:-\d+)?$").match(split):
                    sp = split.split(':')
                    if sp[-1] != "-1":
                        s = slice(int(sp[0]), int(sp[2])+1, int(sp[1])) if len(sp) == 3 else slice(int(sp[0]), int(sp[1])+1)
                    else:
                        s = slice(int(sp[0]), None, int(sp[1])) if len(sp) == 3 else slice(int(sp[0]), None)
                    except_num += tmp[s]
                else:
                    return -1
            for atom in tmp:
                if atom not in except_num:
                    self._in_list_check(atom)
        else:
            return -1
    def _in_list_check(self, atom):
        if atom not in self._atom_list:
            self._atom_list.append(atom)
    def read_atoms_info(self, filepath=None):
        if filepath is None and hasattr(self._args, 'molecules') and path.isfile(self._args.molecules):
            filepath = self._args.molecules
        elif filepath is None and hasattr(self._args, 'molecules') and not path.isfile(self._args.molecules):
            raise FileNotFoundError("'molecules' file not found.")
        elif filepath is None and not hasattr(self._args, 'molecules'):
            raise ValueError("'molecules' argument is not defined.")
        with open(filepath) as read_file:
            tmp = []
            molecule = 0
            for row in read_file:
                row = row.replace('\n', '')
                sp = split(r"[^a-zA-Z0-9]", row)
                if "mol" in sp[0]:
                    tmp.append(sp[1])
                elif compile(r"^\d+(\,\d+)*$").match(row) is not None:
                    for _ in [int(_) for _ in row.split(',')]:
                        self._atoms_info[_]["molecule"] = molecule
                    molecule += 1
    def read_elements(self, software, filepath=None):
        self._elements = software.read_elements() if software.elements == [] else software.elements
        if filepath is not None or (hasattr(self._args, 'sdelements') and self._args.sdelements is not None):
            if filepath is None and path.isfile(self._args.sdelements):
                filepath = self._args.sdelements
            elif filepath is None and not path.isfile(self._args.sdelements):
                raise FileNotFoundError("'sdelements' file not found.")
            with open(filepath) as read_file:
                for idx, line in enumerate(read_file):
                    sp = line.replace('\n', '').split()
                    if idx == 0 and compile(r"^[A-Za-z]+( [A-Za-z]+)*$").match(line) is not None:
                        self._elements = sp
                        break
                    elif compile(r"^[A-Za-z]+ \d+(:\d+){0,2}$").match(line) is not None:
                        elements = [int(num) - self._args.atom for num in sp[1].split(':')]
                        if len(elements) == 1:
                            self._elements[elements[0]] = sp[0]
                        elif len(elements) == 2:
                            for i in range(elements[0], elements[1]+1):
                                self._elements[i] = sp[0]
                        elif len(elements) == 3:
                            for i in range(elements[0], elements[2]+1, elements[1]):
                                self._elements[i] = sp[0]
                    else:
                        print(f"Warning: Input proper format in line {idx+1}. Skip this line.")

def default():
    return 1.9
class Molecule(Fundamental):
    def __init__(self, input_obj=None):
        Fundamental.__init__(self, input_obj)
        del self._AtomStep
        self._mole_dict = defaultdict(list)
        self._mole_kind = []
        self._mole_pos_dict = defaultdict(list)
        self._bond_type = defaultdict(default)
        self._frac_flag, self._cart_flag = False, False
        self._atom_list = None
        if hasattr(input_obj, "bridge"):
            self.bridge = input_obj
        else:
            self.args = input_obj
        if hasattr(self._args, 'molecules') and self._args.molecules is not None:
            self.read_molecules(self._args.molecules)
    @property
    def bridge(self):
        pass
    @property
    def threshold(self):
        return self._bond_type
    @property
    def elements(self):
        return self._elements
    @property
    def lattice(self):
        return self._lattice
    @property
    def molecule_dictionary(self):
        return self._mole_dict
    @property
    def molecule_kind(self):
        return self._mole_kind
    @property
    def molecule_position(self):
        return self._mole_pos_dict
    @property
    def cartesian_flag(self):
        return self._cart_flag
    @property
    def fractional_flag(self):
        return self._frac_flag
    @bridge.setter
    def bridge(self, software):
        self.args, self.lattice, self.elements = software.args, software.lattice, software.elements
        if hasattr(software, 'atom_step'):
            self._atom_list = software.atom_step.atom_list
        elif hasattr(software, 'atom_list'):
            self._atom_list = software.atom_list
    @threshold.setter
    def threshold(self, bond_type):
        if isinstance(bond_type, dict):
            self._bond_type = bond_type
        else:
            raise Exception("'threshold' should be a dictionary.")
    @elements.setter
    def elements(self, elements):
        self._elements = elements
    @lattice.setter
    def lattice(self, lattice):
        self._lattice = lattice
    @molecule_dictionary.setter
    def molecule_dictionary(self, mole_dict):
        if isinstance(mole_dict, dict):
            self._mole_dict = mole_dict
        else:
            raise ValueError("The argument should be 'dict'.")
    @molecule_kind.setter
    def molecule_kind(self, molecule_kind):
        if isinstance(molecule_kind, list):
            self._mole_kind = molecule_kind
        elif isinstance(molecule_kind, dict):
            for key in sorted(list(molecule_kind.keys())):
                self._mole_kind.append(molecule_kind[key])
        else:
            raise ValueError("Molecule kinds should be 'list' or 'dict'.")
    @molecule_position.setter
    def molecule_position(self, MolePosDict):
        if isinstance(MolePosDict, dict):
            self._mole_pos_dict = MolePosDict
        else:
            raise ValueError("The argument should be 'dict'.")
    @cartesian_flag.setter
    def cartesian_flag(self, flag):
        self._cart_flag = flag
    @fractional_flag.setter
    def fractional_flag(self, flag):
        self._frac_flag = flag
    def read_bond(self, filepath=None):
        if self._bond_type == {}:
            if filepath is None and hasattr(self._args, 'bond') and path.isfile(self._args.bond):
                filepath = self._args.bond
            elif filepath is None and hasattr(self._args, 'bond') and not path.isfile(self._args.bond):
                raise FileNotFoundError("'bond' file not found.")
            elif filepath is None and not hasattr(self._args, 'bond'):
                raise ValueError("No 'bond' argument provided.")
            with open(filepath) as read_file:
                row = 0
                for line in read_file:
                    row += 1
                    if len(line.split()[0].split('-')) == 2:
                        if line.split()[0].split('-')[0] not in self._elements:
                            print(f"Warning: {line.split()[0].split('-')[0]} in row {row} does not exist in this system. Skip this row.")
                        elif line.split()[0].split('-')[1] not in self._elements:
                            print(f"Warning: {line.split()[0].split('-')[1]} in row {row} does not exist in this system. Skip this row.")
                        elif compile(r"^\d+\.?\d*$").match(line.split()[1]) is None:
                            print(f"Warning: {line.split()[1]} in row {row} should be a positive number. Skip this row.")
                        else:
                            sp = line.split()[0].split('-')
                            self._bond_type[f'{sp[0]}-{sp[1]}'] = self._bond_type[f'{sp[1]}-{sp[0]}'] = float(line.split()[1])
                    else:
                        print(f"Warning: Input proper format in row {row}. Skip this row.")
    def read_molecules(self, filepath=None):
        if self._mole_dict == {}:
            if filepath is None and hasattr(self._args, 'molecules') and self._args.molecules is not None and path.isfile(self._args.molecules):
                filepath = self._args.molecules
            elif filepath is None and hasattr(self._args, 'molecules') and self._args.molecules is not None and not path.isfile(self._args.molecules):
                raise FileNotFoundError("'molecules' file not found.")
            elif (filepath is None and not hasattr(self._args, 'molecules')) or self._args.molecules is None:
                raise ValueError("No 'molecules' argument provided.")
            with open(filepath) as read_file:
                mole_kind = []
                former_mole_kind = ""
                molecule = 0
                for row in read_file:
                    row = row.replace('\n', '')
                    sp = split(r"[^a-zA-Z0-9]", row)
                    if "mol" in sp[0]:
                        mole_kind.append(sp[1])
                    elif compile(r"^\d+(\,\d+)*$").match(row) is not None:
                        if mole_kind == []:
                            mole_kind.append("mol")
                        if mole_kind[-1] != former_mole_kind:
                            former_mole_kind = mole_kind[-1]
                            kind_counter = 1
                        else:
                            kind_counter += 1
                        if self._atom_list is not None:
                            for atom in sorted([int(_) for _ in row.split(',')]):
                                if atom-self._args.atom in self._atom_list:
                                    self._mole_dict[molecule].append(atom-self._args.atom)
                        else:
                            self._mole_dict[molecule] = sorted([int(_) for _ in row.split(',')])
                        self._mole_kind.append(mole_kind[-1]+f"_{kind_counter}")
                        molecule += 1
    def molecule_wrap(self):
        from kit.accelerate import image_shift, shift_to_origin, distance_matrix_wrap
        if not self._frac_flag:
            from kit.accelerate import c2f_acc
            for key, position in self._mole_pos_dict.items():
                self._mole_pos_dict[key] = c2f_acc(position, self._lattice)
            self._frac_flag = True
        self.build_threshold()
        period_images = array([[0, 0, 1], [1, 0, 1], [-1, 0, 1],
                          [0, 1, 1], [0, -1, 1], [1, 1, 1],
                          [1, -1, 1], [-1, 1, 1], [-1, -1, 1],
                          [0, 0, 0], [1, 0, 0], [-1, 0, 0],
                          [0, 1, 0], [0, -1, 0], [1, 1, 0],
                          [1, -1, 0], [-1, 1, 0], [-1, -1, 0],
                          [0, 0, -1], [1, 0, -1], [-1, 0, -1],
                          [0, 1, -1], [0, -1, -1], [1, 1, -1],
                          [1, -1, -1], [-1, 1, -1], [-1, -1, -1]])
        for key, position in self._mole_pos_dict.items():
            position = shift_to_origin(position)
            first_step_positions = position[0] if position.ndim == 3 else position
            first_step_positions = image_shift(first_step_positions, self._lattice, cartesian=0)
            searched = []
            dist_mat = distance_matrix_wrap(first_step_positions, first_step_positions, self._lattice, period_images)
            for i in range(len(self._mole_dict[key])):
                for j, k in zip(*np.where((dist_mat[i] < max(self._bond_type.values())+0.1) & (dist_mat[i] > 0.01))):
                    if j in searched:
                        continue
                    else:
                        first_step_positions[j] += period_images[k]
                        dist_mat = distance_matrix_wrap(first_step_positions, first_step_positions, self._lattice, period_images)
                        searched.append(j)
                    break
            
            if position.ndim == 3:
                position[0] = first_step_positions
                self._mole_pos_dict[key] = image_shift(position, self._lattice, cartesian=0)
            elif position.ndim == 2:
                position = first_step_positions
    def build_threshold(self):
        from ase.data import atomic_numbers, covalent_radii
        elements = list(set(self._elements))
        for i in range(len(elements)):
            for j in range(i, len(elements)):
                if f"{elements[i]}-{elements[j]}" in self._bond_type.keys() or f"{elements[j]}-{elements[i]}" in self._bond_type.keys():
                    continue
                atomic_num_i, atomic_num_j = atomic_numbers[elements[i]], atomic_numbers[elements[j]]
                bond_threshold = (covalent_radii[atomic_num_i-1] + covalent_radii[atomic_num_j-1]) * 1.2
                self._bond_type[f"{elements[i]}-{elements[j]}"] = self._bond_type[f"{elements[j]}-{elements[i]}"] = bond_threshold if bond_threshold < 2.3 else 2.3

class AtomStep(ABC):
    def __init__(self, input_obj=None, atoms=None, steps=None):
        from kit.accelerate import c2f, f2c
        self._atoms = Atom(input_obj) if atoms is None else atoms
        self._steps = Step(input_obj) if steps is None else steps
        self._atom_list = []
        self._atom_dict = defaultdict(list)
        self._molecules = Molecule(input_obj)
        self._lattice = Lattice()
        self._cart_flag, self._frac_flag, self._fast_flag = False, False, False
        self.c2f, self.f2c = c2f, f2c
        if hasattr(input_obj, "bridge"):
            self.bridge = input_obj
        else:
            self.args = input_obj
    @property
    def args(self):
        return getattr(self, "_args", None)
    @property
    def bridge(self):
        pass
    @property
    def atoms(self):
        return self._atoms
    @property
    def atom_list(self):
        return self._atom_list
    @property
    def elements(self):
        return self._atoms.elements
    @property
    def atoms_info(self):
        return self._atoms.atoms_info
    @property
    def lattice(self):
        return self._lattice.lattice
    @property
    def molecules(self):
        return self._molecules
    @property
    def fast_flag(self):
        return self._fast_flag
    @property
    def fast_position(self):
        if not self._fast_flag:
            raise ValueError("You are not in 'fast' mode.")
        return self._fast_pos
    @property
    def cartesian_flag(self):
        return self._cart_flag
    @property
    def fractional_flag(self):
        return self._frac_flag
    @args.setter
    def args(self, args=None):
        self._args = Args.arg_check(args)
    @bridge.setter
    def bridge(self, software):
        self.args = software.args
        self.atoms = software.atoms
        self._molecules.bridge = software
    @atoms.setter
    def atoms(self, atoms):
        if isinstance(atoms, Atom):
            self._atoms = atoms
            self._atom_list = atoms.get()
        else:
            raise ValueError(f"Only 'Atom' class can be imported (Imported type: {type(atoms)}).")
    @elements.setter
    def elements(self, elements):
        self._atoms.elements = elements
        self._molecules.elements = elements
    @atoms_info.setter
    def atoms_info(self, atoms_info):
        self._atoms.atoms_info = atoms_info
    @lattice.setter
    def lattice(self, lattice):
        lattice = array(lattice)
        if lattice.ndim == 2 or lattice.ndim == 3:
            self._lattice.lattice = lattice
            self._molecules.lattice = lattice
        else:
            raise ValueError("The number of lattice vectors should be 2 or 3.")
    @molecules.setter
    def molecules(self, molecules):
        if isinstance(molecules, Molecule):
            self._molecules = molecules
        else:
            raise ValueError("'molecule' should be 'Molecule' object.")
    @fast_flag.setter
    def fast_flag(self, flag):
        self._fast_flag = flag
        self._fast_pos, self._tmp_pos = [], []
    @fast_position.setter
    def fast_position(self, AtomStep):
        self._fast_pos = AtomStep
    @cartesian_flag.setter
    def cartesian_flag(self, flag):
        self._cart_flag = flag
    @fractional_flag.setter
    def fractional_flag(self, flag):
        self._frac_flag = flag
    def fast_step(self):
        self._fast_pos.append(self._tmp_pos)
        self._tmp_pos = []
    def fast_append(self, atom, AtomStep):
        self._tmp_pos.append(AtomStep)
        self._atoms.elements.append(self._atoms.atoms_info[atom]["element"])
    def put(self, position, atom):
        self._atom_dict[atom].append(position)
    @abstractmethod
    def lock(self):
        if not hasattr(self._args, "molecules") and "molecule" in self._atoms.atoms_info[self._atom_list[0]].keys():
            for key in self._atom_list:
                self._molecules.molecule_dictionary[self._atoms.atoms_info[key]["molecule"]].append(key)
                self._atoms.atoms_info[key].pop("molecule")
                self.build_molecule_position()
        if self._atoms.elements == []:
            for key in self._atom_list:
                self._atoms.elements.append(self._atoms.atoms_info[key]["element"])
                self._atoms.atoms_info[key].pop("element")
        if self._frac_flag:
            self._frac_pos = []
            for atom in self._atoms.get():
                self._frac_pos.append(self._atom_dict[atom])
            self._frac_pos = array(self._frac_pos).transpose((1, 0, 2))
        else:
            self._cart_pos = []
            for atom in self._atoms.get():
                self._cart_pos.append(self._atom_dict[atom])
            self._cart_pos = array(self._cart_pos).transpose((1, 0, 2))
        del self._atom_dict
    def build_molecule_position(self):
        if self._molecules.molecule_dictionary is None or len(self._molecules.molecule_dictionary) == 0:
            self._molecules.read_molecules(self._args.molecules)
        if not self._frac_flag:
            if self._cart_pos.ndim == 3:
                position = self.c2f_acc(self._cart_pos, self.lattice)
            elif self._cart_pos.ndim == 2:
                position = self.c2f(self._cart_pos, self.lattice)
        else:
            position = self._frac_pos
        if position.ndim == 3:
            position = position.transpose((1, 0, 2))
        self._molecules.cartesian_flag, self._molecules.fractional_flag = self._cart_flag, self._frac_flag
        for key, val in self._molecules.molecule_dictionary.items():
            if val[0] not in self._atom_list:
                continue
            for atom in val:
                self._molecules.molecule_position[key].append(position[self._atoms.index_list[atom]])
            self._molecules.molecule_position[key] = array(self._molecules.molecule_position[key]).transpose((1, 0, 2)) if position.ndim == 3 else array(self._molecules.molecule_position[key])
        self._molecules.molecule_wrap()

class AtomStep_Trj(AtomStep):
    def __init__(self, input_obj=None, steps=None, atoms=None):
        AtomStep.__init__(self, input_obj, atoms, steps)
        from kit.accelerate import c2f_acc, f2c_acc
        self.c2f_acc, self.f2c_acc = c2f_acc, f2c_acc
    @property
    def bridge(self):
        pass
    @property
    def steps(self):
        return self._steps
    @property
    def steps_info(self):
        return self._steps.steps_info
    @property
    def cartesian_position(self):
        if self._cart_flag:
            return self._cart_pos
        else:
            self._cart_pos = self.f2c_acc(self._frac_pos, self._lattice.lattice)
            self._cart_flag = True
            return self._cart_pos
    @property
    def fractional_position(self):
        if self._frac_flag:
            return self._frac_pos
        else:
            self._frac_pos = self.c2f_acc(self._cart_pos, self._lattice.lattice)
            self._frac_flag = True
            return self._frac_pos
    @property
    def fast_position(self):
        return self._fast_pos
    @bridge.setter
    def bridge(self, software):
        AtomStep.bridge.fset(self, software)
        self._steps.bridge = software
    @steps_info.setter
    def steps_info(self, steps_info):
        self._steps.steps_info = steps_info
    @cartesian_position.setter
    def cartesian_position(self, cart_pos):
        if cart_pos.ndim == 2:
            from numpy import newaxis
            cart_pos = cart_pos[newaxis, :]
        assert cart_pos.ndim == 3, "Only three-dimensional cartesian data can be imported."
        if not len(cart_pos):
            raise ValueError("The length of 'cartesian_position' cannot be zero.")
        elif isinstance(cart_pos, list):
            self._cart_pos = array(cart_pos)
        elif isinstance(cart_pos, ndarray):
            self._cart_pos = cart_pos
        else:
            raise ValueError("Only 'list' and 'ndarray' can be imported.")
        steps = self._steps.get(slice_flatten=True)
        if self._atom_list == []:
            self._atom_list = self._atoms.get()
        assert self._cart_pos.shape[0] == len(steps), f"The first dimension of cartesian_position should be equal to the number of steps. ({self._cart_pos.shape[0]} != {len(steps)})"
        assert self._cart_pos.shape[1] == len(self._atom_list), f"The second dimension of cartesian_position should be equal to the number of atoms. ({self._cart_pos.shape[1]} != {len(self._atom_list)})"
        self._cart_flag = True
    @fractional_position.setter
    def fractional_position(self, frac_pos):
        if frac_pos.ndim == 2:
            from numpy import newaxis
            frac_pos = frac_pos[:, newaxis]
        assert frac_pos.ndim == 3, "Only three-dimensional fractional data can be imported."
        if not len(frac_pos):
            raise ValueError("The length of 'fractional_position' cannot be zero.")
        elif isinstance(frac_pos, list):
            self._frac_pos = array(frac_pos)
        elif isinstance(frac_pos, ndarray):
            self._frac_pos = frac_pos
        else:
            raise ValueError("Only 'list' and 'ndarray' can be imported.")
        steps = self._steps.get(slice_flatten=True)
        if self._atom_list == []:
            self._atom_list = self._atoms.get()
        assert self._frac_pos.shape[0] == len(steps), f"The first dimension of fractional_position should be equal to the number of steps. ({self._frac_pos.shape[0]} != {len(steps)})"
        assert self._frac_pos.shape[1] == len(self._atom_list), f"The second dimension of fractional_position should be equal to the number of atoms. ({self._frac_pos.shape[1]} != {len(self._atom_list)})"
        self._frac_flag = True
    @fast_position.setter
    def fast_position(self, fast_pos):
        assert fast_pos.ndim == 3, "Only three-dimensional cartesian data can be imported."
        if not len(fast_pos):
            raise ValueError("The length of 'fast_position' cannot be zero.")
        elif isinstance(fast_pos, list):
            self._fast_pos = array(fast_pos)
        elif isinstance(fast_pos, ndarray):
            self._fast_pos = fast_pos
        else:
            raise ValueError("Only 'list' and 'ndarray' can be imported.")
        steps = self._steps.get(slice_flatten=True)
        if self._atom_list == []:
            self._atom_list = self._atoms.get()
        assert self._fast_pos.shape[0] == len(steps), f"The first dimension of fast_position should be equal to the number of steps. ({self._fast_pos.shape[0]} != {len(steps)})"
        assert self._fast_pos.shape[1] == len(self._atom_list), f"The second dimension of fast_position should be equal to the number of atoms. ({self._fast_pos.shape[1]} != {len(self._atom_list)})"
    @steps.setter
    def steps(self, steps):
        if isinstance(steps, Step):
            self._steps = steps
        else:
            raise ValueError("Only 'Step' class can be imported.")
    def lock(self):
        super().lock()

class AtomStep_Single_Point(AtomStep):
    def __init__(self, input_obj=None, atoms=None):
        AtomStep.__init__(self, input_obj, atoms)
        del self._steps
    @property
    def lattice(self):
        return self._lattice.lattice
    @property
    def cartesian_position(self):
        if self._cart_flag:
            return self._cart_pos
        else:
            self._cart_pos = self.f2c(self._frac_pos, self._lattice.lattice)
            self._cart_flag = True
            return self._cart_pos
    @property
    def fractional_position(self):
        if self._frac_flag:
            return self._frac_pos
        else:
            self._frac_pos = self.c2f(self._cart_pos, self._lattice.lattice)
            self._frac_flag = True
            return self._frac_pos
    @property
    def fast_position(self):
        return self._fast_pos
    @lattice.setter
    def lattice(self, lattice):
        lattice = array(lattice)
        if lattice.ndim == 2:
            self._lattice.lattice = lattice
            self._molecules.lattice = lattice
        elif lattice.ndim == 3:
            self._lattice.lattice = lattice[0]
            self._molecules.lattice = lattice[0]
        else:
            raise ValueError("The number of lattice vectors should be 2.")
    @cartesian_position.setter
    def cartesian_position(self, cart_pos):
        if isinstance(cart_pos, list) and not len(cart_pos):
            self._cart_pos = array(cart_pos)
        elif isinstance(cart_pos, list):
            self._cart_pos = array(cart_pos)
        elif isinstance(cart_pos, ndarray):
            self._cart_pos = cart_pos
        else:
            raise ValueError("Only 'list' and 'ndarray' can be imported.")
        assert self._cart_pos.ndim == 2, "Only two-dimensional cartesian data can be imported."
        if self._atom_list == []:
            self._atom_list = self._atoms.get()
        assert self._cart_pos.shape[0] == len(self._atom_list), f"The first dimension of cartesian_position should be equal to the number of atoms. ({self._cart_pos.shape[0]} != {len(self._atom_list)})"
        self._cart_flag = True
    @fractional_position.setter
    def fractional_position(self, frac_pos):
        if isinstance(frac_pos, list) and not len(frac_pos):
            self._frac_pos = array(frac_pos)
        elif isinstance(frac_pos, list):
            self._frac_pos = array(frac_pos)
        elif isinstance(frac_pos, ndarray):
            self._frac_pos = frac_pos
        else:
            raise ValueError("Only 'list' and 'ndarray' can be imported.")
        assert self._frac_pos.ndim == 2, f"Only two-dimensional fractional data can be imported. ({frac_pos.ndim} != 2)"
        if self._atom_list == []:
            self._atom_list = self._atoms.get()
        assert self._frac_pos.shape[0] == len(self._atom_list), f"The first dimension of fractional_position should be equal to the number of atoms. ({self._frac_pos.shape[0]} != {len(self._atom_list)})"
        self._frac_flag = True
    @fast_position.setter
    def fast_position(self, fast_pos):
        if isinstance(fast_pos, list) and not len(fast_pos):
            self._fast_pos = array(fast_pos)
        elif isinstance(fast_pos, list):
            self._fast_pos = array(fast_pos)
        elif isinstance(fast_pos, ndarray):
            self._fast_pos = fast_pos
        else:
            raise ValueError("Only 'list' and 'ndarray' can be imported.")
        assert self._fast_pos.ndim == 2, "Only two-dimensional fractional data can be imported."
        if self._atom_list == []:
            self._atom_list = self._atoms.get()
        assert self._fast_pos.shape[0] == len(self._atom_list), f"The first dimension of fast_position should be equal to the number of atoms. ({self._fast_pos.shape[0]} != {len(self._atom_list)})"
    def lock(self):
        super().lock()

class Convert:
    @staticmethod
    def elements2eleNum(elements):
        eleNum = {}
        elementOrder = []
        for element in elements:
            if element not in eleNum.keys():
                eleNum[element] = 1
                elementOrder.append(element)
            else:
                eleNum[element] += 1
        return eleNum, elementOrder
    @staticmethod
    def eleNum2elements(ele, num):
        elements = []
        for idx, element in enumerate(ele):
            for _ in range(int(num[idx])):
                elements.append(element)
        return elements

class Args:
    def __init__(self, default=None, arg_func=None):
        args = ArgumentParser()
        if arg_func is None:
            arg_func = self.func
        
        args = arg_func(default, args)
        
        self.__args = args.parse_args()
    def func(self, default=None, args=None):
        args.add_argument("input", type=str, nargs='?', help="input file or directory path")
        args.add_argument("-step", dest="step", type=int, default=1, help="initial step number of the trajectory file. default is 1")
        args.add_argument("-atom", dest="atom", type=int, default=0, help="initial atom number. default is 0")
        args.add_argument("-stepfile", dest="stepfile", type=str, default=None, help="provide the step information. default is None")
        args.add_argument("-atomfile", dest="atomfile", type=str, default=None, help="provide the atom information. default is None")
        args.add_argument("-elementfile", dest="elementfile", type=str, default="Atoms.txt", help="provide the elements information when the input file does not include, e.g., QE, and CONQUEST. default is Atoms.txt.")
        args.add_argument("-sdelements", dest="sdelements", type=str, default=None, help="provide self-defined elements information.")
        if isinstance(default, dict):
            for key in default.keys():
                if key == "output" and not isinstance(default[key], list):
                    args.add_argument("-output", dest="output", type=str, default=default[key], help=f"output file or directory name. default is {default[key]}")
                elif isinstance(default[key], Iterable) and len(default[key]) == 3:
                    args.add_argument(f"-{key}", dest=key, type=default[key][1], default=default[key][0], help=f"{default[key][2]}")
                elif isinstance(default[key], Iterable) and len(default[key]) == 2:
                    args.add_argument(f"-{key}", dest=key, type=default[key][0].__class__, default=default[key][0], help=f"{default[key][1]}")
                elif isinstance(default[key], Iterable) and isinstance(default[key], str):
                    args.add_argument(f"-{key}", dest=key, type=default[key].__class__, default=default[key], help=f"{default[key]}")
                elif isinstance(default[key], Iterable):
                    args.add_argument(f"-{key}", dest=key, type=default[key][0].__class__, default=default[key][0], help=f"{default[key][0]}")
                else:
                    args.add_argument(f"-{key}", dest=key, type=default[key].__class__, default=default[key], help=f"{key}")
        elif default is not None:
            raise ValueError("default should be a dictionary.")
        return args
    @property
    def args(self):
        return self.__args
    @property
    def input_type(self):
        return self.__inputType
    @input_type.setter
    def input_type(self, inputType):
        self.__inputType = inputType
    @staticmethod
    def same_name(present_path, filename, factor=0):
        while(1):
            if "_Alpha=" in filename:
                repeat = True
            else:
                repeat = False
            if filename is None:
                break
            if factor and not repeat:
                tmp = path.splitext(filename)
                filename = tmp[0] + f"_Alpha={factor}" + tmp[1]
            file = path.join(present_path, filename)
            if path.isfile(file):
                action = input(f"Warning: '{filename}' file exists in '{present_path}'. Delete, move, or rename (d/m/n): ")
                repeat = True
            elif path.isfile(f"{file}.csv"):
                action = input(f"Warning: '{filename}.csv' file exists in '{present_path}'. Delete, move, or rename (d/m/n): ")
                filename += ".csv"
                repeat = True
            elif path.isfile(f"{file}.dat"):
                action = input(f"Warning: '{filename}.dat' file exists in '{present_path}'. Delete, move, or rename (d/m/n): ")
                filename += ".dat"
                repeat = True
            elif path.isdir(file):
                action = input(f"Warning: '{filename}' directory exists in '{present_path}'. Delete, move, or rename (d/m/n): ")
                repeat = True
            else:
                return filename
            Args.action(present_path, filename, action)
    def input_check(self):
        self.__args.input = self.input_file_check(self.__args.input, self.__inputType)
    def action(present_path, filename, action):
        filePath = path.join(present_path, filename)
        if action.lower() == 'd' and path.isdir(filePath):
            from shutil import rmtree
            rmtree(filePath)
        elif action.lower() == 'd' and path.isfile(filePath):
            from os import remove
            remove(filePath)
        elif action.lower() == 'n':
            from os import rename
            newName = input("New name: ")
            rename(filePath, path.join(present_path, newName))
        elif action.lower() == 'm':
            from shutil import move
            newPath = input("New path: ")
            if path.isdir(newPath):
                move(filePath, newPath)
            else:
                print("Warning: Provide a path to existence.")
    @staticmethod
    def input_file_check(input_path, input_type, mandatory=True, isdir=False):
        while(1):
            if input_path is None:
                input_path = input(f"Input {input_type} path: ")
            elif input_path is None and not mandatory:
                input_path = input(f"Input {input_type} path (If not needed, input 'no'): ")
            if not mandatory and input_path.lower() == "no":
                return None
            if not isdir:
                if path.isfile(input_path):
                    return input_path
                elif path.isfile(path.join(input_path, input_type)):
                    return path.join(input_path, input_type)
                elif path.isfile(input_path+input_type):
                    return input_path+input_type
                elif path.isfile(input_path+"."+input_type):
                    return input_path+"."+input_type
                else:
                    print(f"Warning: '{input_type}' does not exist in {path.abspath(input_path)}.\n")
                    input_path = None
            else:
                if path.isdir(input_path):
                    return input_path
                elif path.isdir(path.join(input_path, input_type)):
                    return path.join(input_path, input_type)
                else:
                    print(f"Warning: {input_type} does not exist in {path.abspath(input_path)}.")
    @staticmethod
    def arg_check(args):
        if isinstance(args, dict):
            args = Args(args)
        if isinstance(args, Args):
            return args.args
        elif isinstance(args, Namespace):
            return args
        else:
            return None

class Lattice:
    def __init__(self):
        self.__first_flag, self.__NpT_flag = False, False
        self.__firstLattice = []
    @property
    def NpT_flag(self):
        return self.__NpT_flag
    @property
    def lattice(self):
        if isinstance(self.__firstLattice, list):
            return self.__firstLattice
        if self._lattice is None:
            return array(self.__firstLattice)
        return array(self._lattice)
    @NpT_flag.setter
    def NpT_flag(self, NpT_flag):
        self.__NpT_flag = NpT_flag
    @lattice.setter
    def lattice(self, lattice):
        if not self.__first_flag:
            self.__first_flag = True
            self.__firstLattice = array(lattice)
            self._lattice = None
        elif self._lattice is None and np.allclose(self.__firstLattice, array(lattice)):
            self._lattice = lattice
        elif self._lattice is None:
            self.__NpT_flag = True
            self._lattice = [self.__firstLattice]
            self._lattice.append(lattice)
        elif self.__NpT_flag:
            self._lattice.append(array(lattice))

class Graph(Position, Periodic):
    def __init__(self, input_obj=None):
        Position.__init__(self, input_obj)
        Periodic.__init__(self, input_obj)
        self._dist_mat = None
        self._adj_mat = None
        self._UnupdatedBondList = []
        self._period_images = array([[0, 0, 1], [1, 0, 1], [-1, 0, 1],
                                     [0, 1, 1], [0, -1, 1], [1, 1, 1],
                                     [1, -1, 1], [-1, 1, 1], [-1, -1, 1],
                                     [0, 0, 0], [1, 0, 0], [-1, 0, 0],
                                     [0, 1, 0], [0, -1, 0], [1, 1, 0],
                                     [1, -1, 0], [-1, 1, 0], [-1, -1, 0],
                                     [0, 0, -1], [1, 0, -1], [-1, 0, -1],
                                     [0, 1, -1], [0, -1, -1], [1, 1, -1],
                                     [1, -1, -1], [-1, 1, -1], [-1, -1, -1]])
        if hasattr(input_obj, "bridge"):
            self.bridge = input_obj
        from kit.accelerate import distance_matrix_wrap, distance_matrix
        self.distance_matrix_wrap, self.distance_matrix_func = distance_matrix_wrap, distance_matrix
    @property
    def atom_step(self):
        return self._AtomStep
    @property
    def bridge(self):
        pass
    @property
    def period_images(self):
        return self._period_images
    @property
    def distance_matrix(self):
        return self._dist_mat
    @property
    def adjacent_matrix(self):
        return self._adj_mat
    @property
    def adjacent_list(self):
        return self._AdjList
    @property
    def molecules(self):
        return self._AtomStep.molecules
    @atom_step.setter
    def atom_step(self, atom_step):
        if isinstance(atom_step, AtomStep_Trj) or isinstance(atom_step, AtomStep_Single_Point):
            self._AtomStep = atom_step
            self.molecules = self._AtomStep.molecules
        else:
            raise ValueError("'atom_step' should be ' 'AtomStep_Trj' or 'AtomStep_Single_Point'.")
    @bridge.setter
    def bridge(self, software):
        self.args, self.atom_step = software.args, software.atom_step
        if hasattr(self._args, "molecules") and self._AtomStep.molecules.molecule_dictionary != {}:
            self.read_molecules(self._args.molecules)
    @period_images.setter
    def period_images(self, period_images):
        self._period_images = array(period_images)
        self._period_images_length = len(self._period_images)
    @distance_matrix.setter
    def distance_matrix(self, dist_mat):
        if isinstance(dist_mat, ndarray):
            self._dist_mat = dist_mat
        else:
            raise ValueError("'distance_matrix' should be a numpy array.")
    @adjacent_matrix.setter
    def adjacent_matrix(self, adj_mat):
        self._adj_mat = adj_mat
    @adjacent_list.setter
    def adjacent_list(self, adjList):
        if isinstance(self._AdjList, list):
            self._AdjList = adjList
        else:
            raise ValueError("'adjacent_list' should be a dictionary.")
    @molecules.setter
    def molecules(self, molecules):
        self._AtomStep.molecules = molecules
        self._fragment_mole_dict = deepcopy(self._AtomStep.molecules.molecule_dictionary)
        self._fragment_mole_pos = deepcopy(self._AtomStep.molecules.molecule_position)
    def graph_to_molecule_position(self, graph=None, molecule_dictionary=None, molecule_position=None, wrap=True):
        graph = self._adj_mat if graph is None and isinstance(self._adj_mat, ndarray) else graph
        mole_dict = self._AtomStep.molecules.molecule_dictionary if molecule_dictionary is None and isinstance(self._AtomStep.molecules.molecule_dictionary, dict) else molecule_dictionary
        mole_pos = self._AtomStep.molecules.molecule_position if molecule_position is None and isinstance(self._AtomStep.molecules.molecule_position, dict) else molecule_position
        atoms = self._AtomStep.atom_list
        mol_num = 0
        for i, j in zip(*np.where(graph > 0)):
            new_id, searched_id = [], []
            for molecule in mole_dict.values():
                if atoms[j] in molecule:
                    break
            else:
                for neighbor_atom in np.where(graph[j] > 0)[0]:
                    if neighbor_atom not in new_id and neighbor_atom not in searched_id:
                        new_id.append(neighbor_atom)
                if new_id == []:
                    continue
                while(1):
                    poped = new_id.pop(0)
                    searched_id.append(poped)
                    for neighbor_atom in np.where(graph[poped])[0]:
                        if neighbor_atom not in new_id and neighbor_atom not in searched_id:
                            new_id.append(neighbor_atom)
                    if new_id == []:
                        mole_dict[mol_num].extend([atoms[_] for _ in searched_id])
                        break
                mol_num += 1
        
        for atom in atoms:
            for mole_atoms in mole_dict.values():
                if atom in mole_atoms:
                    break
            else:
                mole_dict[mol_num].append(atom)
                mol_num += 1
        positions = self._AtomStep.fractional_position
        for mol_num, molecule in mole_dict.items():
            molecule.sort()
            atoms_idx = self._AtomStep.atoms.index_list
            if isinstance(mole_pos[mol_num], list) and mole_pos[mol_num] == []:
                if positions.ndim == 2:
                    for atom in molecule:
                        mole_pos[mol_num].append(positions[atoms_idx[atom]])
                    mole_pos[mol_num] = array(mole_pos[mol_num])
                elif positions.ndim == 3:
                    for atom in molecule:
                        mole_pos[mol_num].append(positions[:, atoms_idx[atom]])
                    mole_pos[mol_num] = array(mole_pos[mol_num]).transpose((1, 0, 2))
        if wrap:
            self._AtomStep.molecules.fractional_flag = True
            self._AtomStep.molecules.molecule_wrap()
        return mole_dict, mole_pos
    def build_graph(self, wrap=True, molecule=True):
        elements = self._AtomStep.atoms.elements
        self._adj_mat = np.zeros((len(elements), len(elements)))
        if molecule:
            if self._AtomStep.molecules.molecule_position == {}:
                self._AtomStep.build_molecule_position()
            self.molecule_distance_matrix()
        else:
            positions = self._AtomStep.fractional_position if self._AtomStep.fractional_position.ndim == 2 else self._AtomStep.fractional_position[0]
            if positions.shape[0] < 1000:
                self._dist_mat = self.build_distance_matrix(wrap)
            else:
                self._dist_mat = np.zeros((positions.shape[0], positions.shape[0]))
                for idx, position in enumerate(positions):
                    self._dist_mat[idx] = self.distance_matrix_wrap(position[np.newaxis, :], positions, self._AtomStep.lattice, self._wrap).min(axis=2) if wrap else self.distance_matrix_func(position[np.newaxis, :], positions, self._AtomStep.lattice)
            self._adj_mat, self._AdjList = self.build_graph_from_distance_matrix(self._dist_mat)
    def build_distance_matrix(self, wrap=True):
        position = self._AtomStep.fractional_position if self._AtomStep.fractional_position.ndim == 2 else self._AtomStep.fractional_position[0]
        if wrap:
            return self.distance_matrix_wrap(position, position, self._AtomStep.lattice, self._period_images)
        else:
            return self.distance_matrix_func(position, position, self._AtomStep.lattice)
    def molecule_distance_matrix(self):
        atoms_idx = self._AtomStep.atoms.index_list
        elements = self._AtomStep.atoms.elements
        for mole_num, MolePos in self._AtomStep.molecules.molecule_position.items():
            if len(MolePos) == 1:
                continue
            mole_dict = self._AtomStep.molecules.molecule_dictionary[mole_num]
            position = MolePos if MolePos.ndim == 2 else MolePos[0]
            dist_mat = self.distance_matrix_func(position, position, self._AtomStep.lattice)
            for i, j in zip(*np.where((dist_mat < max(self._AtomStep.molecules.threshold.values())) & (dist_mat > 0.01))):
                atom_idx_i, atom_idx_j = atoms_idx[mole_dict[i]], atoms_idx[mole_dict[j]]
                bond_type = f"{elements[atom_idx_i]}-{elements[atom_idx_j]}"
                threshold = self._AtomStep.molecules.threshold[bond_type]
                if dist_mat[i, j] <= threshold:
                    self._adj_mat[atoms_idx[mole_dict[i]]][atoms_idx[mole_dict[j]]] = self._adj_mat[atoms_idx[mole_dict[j]]][atoms_idx[mole_dict[i]]] = threshold
            atom_idx_i, atom_idx_j = np.unravel_index(np.argmax(dist_mat), dist_mat.shape)
            while(len(self.connect(mole_dict[atom_idx_i], mole_dict[atom_idx_j], self._adj_mat)) > 1 and len(mole_dict) > 1):
                fragment_i = self.connect(mole_dict[atom_idx_i], mole_dict[atom_idx_j], self._adj_mat)[0]
                fragment_j = [atom for atom in mole_dict if atom not in fragment_i]
                fragment_pos_i, fragment_pos_j = [], []
                for i in fragment_i:
                    fragment_idx = mole_dict.index(i)
                    fragment_pos_i.append(position[fragment_idx])
                for j in fragment_j:
                    fragment_idx = mole_dict.index(j)
                    fragment_pos_j.append(position[fragment_idx])
                dist_mat = self.distance_matrix_func(array(fragment_pos_i), array(fragment_pos_j), self._AtomStep.lattice)
                fragment_atom_idx_i, fragment_atom_idx_j = np.unravel_index(np.argmin(dist_mat), dist_mat.shape)
                bond_type = f"{elements[atoms_idx[mole_dict[fragment_atom_idx_i]]]}-{elements[atoms_idx[mole_dict[fragment_atom_idx_j]]]}"
                self._adj_mat[atoms_idx[fragment_j[fragment_atom_idx_j]], atoms_idx[fragment_i[fragment_atom_idx_i]]] = self._adj_mat[atoms_idx[fragment_i[fragment_atom_idx_i]], atoms_idx[fragment_j[fragment_atom_idx_j]]] = self._AtomStep.molecules.threshold[bond_type]
    def build_graph_from_distance_matrix(self, distance_matrix=None):
        dist_mat = distance_matrix if distance_matrix is not None and isinstance(distance_matrix, ndarray) else self._dist_mat
        dist_mat = array(dist_mat).min(axis=2) if dist_mat.ndim == 3 else dist_mat
        self._AtomStep.molecules.build_threshold()
        threshold = max(self._AtomStep.molecules.threshold.values())
        adj_mat = np.zeros(dist_mat.shape)
        elements = self._AtomStep.atoms.elements
        if dist_mat.shape[0] > 3000:
            index_iter = ((i, j) for i in range(adj_mat.shape[0]) for j in np.where((dist_mat[i] < threshold) & (dist_mat[i] > 0.01))[0])
        else:
            index_iter = zip(*np.where((dist_mat < threshold) & (dist_mat > 0.01)))
        for i, j in index_iter:
            if adj_mat[i, j] > 0:
                pass
            else:
                bond_type = f"{elements[i]}-{elements[j]}"
                adj_mat[i, j] = adj_mat[j, i] = self._AtomStep.molecules.threshold[bond_type] if dist_mat[i, j] < self._AtomStep.molecules.threshold[bond_type] else 0
        return adj_mat, self.matrix_to_list(adj_mat)
    def build_graph_from_molecule_dictionary(self, molecule_dictionary=None, distance_matrix=None):
        mole_dict = molecule_dictionary if molecule_dictionary is not None and isinstance(molecule_dictionary, dict) else self._AtomStep.molecules.molecule_dictionary
        dist_mat = distance_matrix if distance_matrix is not None and isinstance(distance_matrix, ndarray) else self._dist_mat
        elements = self._AtomStep.atoms.elements
        adj_mat = np.zeros((dist_mat.shape[0], dist_mat.shape[0]))
        atoms_idx = self._AtomStep.atoms.index_list
        for molecule in mole_dict.values():
            mole_idx = [atoms_idx[atom] for atom in molecule]
            for i in mole_idx:
                dist_arr = [dist_mat[i][idx] for idx in mole_idx]
                j = dist_mat[i].index(min(dist_arr[dist_arr > 0]))
                bond_type = f"{elements[i]}-{elements[j]}"
                adj_mat[i, j] = adj_mat[j, i] = self._AtomStep.molecules.threshold[bond_type]
        return adj_mat
    def matrix_to_list(self, adjacent_matrix=None):
        adj_mat = adjacent_matrix if adjacent_matrix is not None and isinstance(adjacent_matrix, ndarray) else self._adj_mat
        AdjList = defaultdict(list)
        for i, j in zip(*np.where(adj_mat > 0)):
            atom_i, atom_j = self._AtomStep.atom_list[i], self._AtomStep.atom_list[j]
            AdjList[atom_i].append([atom_j, adj_mat[i, j]])
            AdjList[atom_j].append([atom_i, adj_mat[j, i]])
        return AdjList
    def list_to_matrix(self, adjacent_list=None):
        AdjList = adjacent_list if adjacent_list is not None and isinstance(adjacent_list, dict) else self._AdjList
        adj_mat = np.zeros((len(self._AtomStep.atoms.get()), len(self._AtomStep.atoms.get())))
        atoms_idx = self._AtomStep.atoms.index_list
        for cen_atom, mea_atoms in AdjList.items():
            cen_atom_idx = atoms_idx[cen_atom]
            for mea_atom, bond_length in mea_atoms:
                mea_atom_idx = atoms_idx[mea_atom]
                adj_mat[cen_atom_idx, mea_atom_idx] = adj_mat[mea_atom_idx, cen_atom_idx] = bond_length
        return adj_mat
    def create_bond(self, atom1, atom2, graph=None):
        elements = self._AtomStep.atoms.elements
        graph = self._adj_mat if graph is None else graph
        atom1_idx, atom2_idx = self._AtomStep.atoms.index_list[atom1], self._AtomStep.atoms.index_list[atom2]
        element1, element2 = elements[atom1_idx], elements[atom2_idx]
        graph[atom1_idx, atom2_idx] = graph[atom2_idx, atom1_idx] = self._AtomStep.molecules.threshold[f"{element1}-{element2}"]
        if [atom2, atom1, 1] not in self._UnupdatedBondList:
            self._UnupdatedBondList.append([atom1, atom2, 1])
        return graph
    def delete_bond(self, atom1, atom2, graph=None):
        graph = self._adj_mat if graph is None else graph
        atom1_idx, atom2_idx = self._AtomStep.atoms.index_list[atom1], self._AtomStep.atoms.index_list[atom2]
        graph[atom1_idx, atom2_idx] = graph[atom2_idx, atom1_idx] = 0
        if [atom2, atom1, 0] not in self._UnupdatedBondList:
            self._UnupdatedBondList.append([atom1, atom2, 0])
        return graph
    def update_graph(self):
        mole_dict, mole_pos_dict = self._fragment_mole_dict, self._fragment_mole_pos
        while(self._UnupdatedBondList != []):
            bond = self._UnupdatedBondList.pop(0)
            for mol, mol_atoms in mole_dict.items():
                if bond[0] in mol_atoms or bond[1] in mol_atoms:
                    molecule = mol_atoms
                    break
            fragments = self.connect(bond[0], bond[1], self._adj_mat, molecule)
            new_key = max(mole_dict.keys())+1
            if bond[2]:
                mole_num = [k for k, v in mole_dict.items() if bond[0] in v or bond[1] in v]
                mole_dict[new_key] = mole_dict[mole_num[0]] + mole_dict[mole_num[1]]
                del self._fragment_mole_dict[mole_num[0]], self._fragment_mole_dict[mole_num[1]]
                pos0 = mole_pos_dict[mole_num[0]].transpose((1, 0, 2))
                pos1 = mole_pos_dict[mole_num[1]].transpose((1, 0, 2))
                mole_pos_dict[new_key] = np.concatenate([pos0, pos1], axis=0).transpose((1, 0, 2))
                del self._fragment_mole_pos[mole_num[0]], self._fragment_mole_pos[mole_num[1]]
            elif len(fragments) > 1:
                positions = mole_pos_dict[mol].transpose((1, 0, 2))
                atom_to_pos = dict(zip(mol_atoms, positions))
                for i, fragment in enumerate(fragments):
                    if fragment in mole_dict.values():
                        continue
                    atoms_in_frag = [atom for atom in fragment if atom in atom_to_pos]
                    mole_dict[new_key + i] = atoms_in_frag
                    mole_pos_dict[new_key + i] = array([atom_to_pos[atom] for atom in atoms_in_frag]).transpose((1, 0, 2))
                del self._fragment_mole_dict[mol], self._fragment_mole_pos[mol]
    def connect(self, atom1, atom2, graph=None, molecule=None):
        graph = self._adj_mat if graph is None else graph
        atom1_idx, atom2_idx = self._AtomStep.atoms.index_list[atom1], self._AtomStep.atoms.index_list[atom2]
        searched = []
        mole_list = self.walking(atom1_idx, graph)
        searched.append(sorted([self._AtomStep.atom_list[atom] for atom in mole_list]))
        if atom2 not in searched[0]:
            mole_list = self.walking(atom2_idx, graph)
            searched.append(sorted([self._AtomStep.atom_list[atom] for atom in mole_list]))
        if molecule is not None:
            tmp = [self._AtomStep.atoms.index_list[atom] for atom in molecule if atom not in searched[0] and atom not in searched[1]]
            while(tmp != []):
                atom = tmp.pop(0)
                mole_list = self.walking(atom, graph)
                searched.append(sorted([self._AtomStep.atom_list[atom] for atom in mole_list]))
                for atom in mole_list:
                    if atom in tmp:
                        tmp.remove(atom)
        return searched
    @staticmethod
    def walking(index, graph):
        new, mole_list = [index], []
        while(new != []):
            atom = new.pop(0)
            if atom not in mole_list:
                mole_list.append(atom)
            connected_atom = list(np.where(graph[atom])[0])
            new.extend([atom for atom in connected_atom if atom not in mole_list])
        return mole_list
    @staticmethod
    def walking_list(index, graph):
        new, mole_list = [index], []
        while(new != []):
            atom = new.pop(0)
            if atom not in mole_list:
                mole_list.append(atom)
            connected_atom = list(atom2 for atom2 in graph[atom].keys() if graph[atom][atom2])
            new.extend([atom for atom in connected_atom if atom not in mole_list])
        return mole_list
