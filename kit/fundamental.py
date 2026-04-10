from abc import ABC, abstractmethod
from collections.abc import Iterable
from os import path
from copy import deepcopy
from argparse import Namespace
from collections import defaultdict
from re import fullmatch, split
from numpy import ndarray, array
import numpy as np
from kit.args import Args

class Fundamental(ABC):
    def __init__(self, input_obj=None, steps=None, atoms=None, **kwargs):
        super().__init__(**kwargs)
        self._title = "title"
        self._AtomStep = self._create_atom_step(input_obj, steps, atoms, **kwargs)
        self._elements = []
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
        return self._elements
    @property
    def atoms_info(self):
        return self._AtomStep.atoms.atoms_info
    @property
    def atom_step(self):
        return self._AtomStep
    @property
    def Lattice(self):
        return self._AtomStep.Lattice
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
        self._args = Args.arg_check(args)
        if hasattr(self._args, "input"):
            self._title = self._args.input.split('.')[0] if hasattr(self._args, 'input') and self._args.input is not None else "position"
    @title.setter
    def title(self, title):
        self._title = str(title)
    @elements.setter
    def elements(self, elements):
        self._elements = elements
    @atoms_info.setter
    def atoms_info(self, atoms_info):
        self._AtomStep.atoms.atoms_info = atoms_info
    @atom_step.setter
    def atom_step(self, atom_step):
        if atom_step.__class__.__bases__[0] == AtomStep:
            self._AtomStep = atom_step
        else:
            raise ValueError("Only 'AtomStep_Trj' or 'AtomStep_Single_Point' can be imported.")
    @Lattice.setter
    def Lattice(self, Lattice):
        self._AtomStep.Lattice = Lattice
    @lattice.setter
    def lattice(self, lattice):
        self._AtomStep.lattice = lattice
    @atoms.setter
    def atoms(self, atoms):
        self._AtomStep.atoms = atoms
    @molecules.setter
    def molecules(self, molecules):
        self._AtomStep.molecules = molecules
    def read_molecules(self, filepath=None):
        if filepath is not None and not path.isfile(filepath):
            filepath = None
        if hasattr(self._args, "molecules") and self._args.molecules is not None and not path.isfile(self._args.molecules):
            self._args.molecules = None
        if filepath is None:
            if not hasattr(self._args, "molecules"):
                raise ValueError("'molecules' argument does not exist.")
            elif self._args.molecules is None:
                while(1):
                    filepath = input("Input the molecules file path: ")
                    if path.isfile(filepath):
                        with open(filepath) as read_file:
                            for row in read_file:
                                row = row.replace('\n', '')
                                if fullmatch(r"mol [a-zA-Z0-9]+", row) is None and fullmatch(r"\d+(\,\d+)*", row) is None:
                                    print("Warning: Format of the molecules file is wrong.")
                                    break
                            else:
                                self._args.molecules = filepath
                                break
                    else:
                        print("Warning: File not found.")
            else:
                filepath = self._args.molecules
        
        self._AtomStep.molecules.read_molecules(filepath)
        
        self.atoms.molecule_kind = self._AtomStep.molecules.molecule_kind
        self.atoms.molecule_dictionary = self._AtomStep.molecules.molecule_dictionary
        for atom_idx, molecule in self._AtomStep.molecules.molecule_dictionary.items():
            for atom in molecule:
                self._AtomStep.atoms.atoms_info[atom]["molecule"] = atom_idx
    def _create_atom_step(self, input_obj=None, steps=None, atoms=None, **kwargs):
        return None

class Position(Fundamental):
    def __init__(self, input_obj=None, **kwargs):
        super().__init__(input_obj=input_obj, **kwargs)
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
        if fullmatch(r"\d+(\.\d+)?", str(scaling)) is not None:
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
    def __init__(self, input_obj=None, **kwargs):
        super().__init__(input_obj=input_obj, **kwargs)
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
    def __init__(self, input_obj=None, steps=None, atoms=None, **kwargs):
        self._step_counter = 0
        self._cutoff_step = -1
        super().__init__(input_obj=input_obj, steps=steps, atoms=atoms, **kwargs)
        self._title = "trajectory"
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
        return self._step_counter + self._args.step - 1
    @Fundamental.atom_step.setter
    def atom_step(self, AtomStep):
        if isinstance(AtomStep, AtomStep_Trj):
            self._AtomStep = deepcopy(AtomStep)
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
    def _create_atom_step(self, input_obj=None, steps=None, atoms=None, **kwargs):
        return AtomStep_Trj(input_obj, steps, atoms)

class Single_Point(Periodic, Position):
    def __init__(self, input_obj=None, atoms=None, **kwargs):
        self._relax_flag = False
        self._relax_border = []
        super().__init__(input_obj=input_obj, atoms=atoms, **kwargs)
        self._AtomStep = self._create_atom_step(input_obj=input_obj, atoms=atoms, **kwargs)
    @property
    def relax(self):
        return self._relax_flag
    @property
    def relax_border(self):
        return self._relax_border
    @Fundamental.atom_step.setter
    def atom_step(self, AtomStep):
        if isinstance(AtomStep, AtomStep_Single_Point):
            self._AtomStep = deepcopy(AtomStep)
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
        if isinstance(relax_border, Iterable):
            self._relax_border = relax_border
        elif fullmatch(r"-?\d+(\.\d+)?", str(relax_border)) is not None:
            relax_border = float(relax_border)
            if -1 <= relax_border <= 1:
                self._relax_border = relax_border
            else:
                raise ValueError("The relax border should be between -1 and 1")
        else:
            raise ValueError("The relax border should give a list or a number between -1 and 1")
    def _create_atom_step(self, input_obj=None, steps=None, atoms=None, **kwargs):
        return AtomStep_Single_Point(input_obj, atoms)

class Charge:
    def __init__(self, **kwargs):
        self._charge, self._diff, self._ref = [], [], {}
        super().__init__(**kwargs)
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

def _none_factory():
    return None

class Step:
    def __init__(self, input_obj=None, put_flag=False, **kwargs):
        super().__init__(**kwargs)
        self.__steps = []
        self._steps_info = defaultdict(dict)
        self._total_steps = -1
        self._index_dict = defaultdict(_none_factory)
        if hasattr(input_obj, "bridge"):
            self.bridge = input_obj
            if put_flag:
                if isinstance(input_obj, Step):
                    self.put(input_obj)
                else:
                    self.put(input_obj.steps)
        elif isinstance(input_obj, Args) or isinstance(input_obj, Namespace) or isinstance(input_obj, dict):
            self._args = Args.arg_check(input_obj)
            self._start = self._args.step if hasattr(self._args, "step") else 1
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
            for index, step in enumerate(self.__steps):
                self._index_dict[step] = index
        return self._index_dict
    @start.setter
    def start(self, start):
        if fullmatch(r"\-?\d+", str(start)) is not None:
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
            for step, index in index_list.items():
                self._index_dict[step] = index
        else:
            raise ValueError("'index_list' should be 'dict'.")
    def put(self, input_steps):
        if fullmatch(r"\-?\d+", str(input_steps)) is not None:
            input_steps = int(input_steps)
            if input_steps > self._total_steps:
                return -2
            elif input_steps < 0 and input_steps+self._total_steps+1 < 0:
                return -3
            self._steps(input_steps)
        elif ':' in input_steps:
            if fullmatch(r"\-?\d+:\-?\d+(:\-?\d+)?", str(input_steps)) is None:
                return -4
            sp = [int(step) for step in input_steps.split(':')]
            if (len(sp) == 3 and int(sp[0]) > int(sp[2])) or (len(sp) == 2 and int(sp[0]) > int(sp[1])):
                return -5
            elif min(sp) < 0 and min(sp)+self._total_steps+1 < 0:
                return -3
            elif max(sp) > self._total_steps:
                return -2
            else:
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
        elif str(input_steps).isalpha():
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
        self.__steps = sorted(self.__steps)
    @property
    def err_list(self):
        return defaultdict(int, {-5: "Warning: The start number should be less than the end number.", \
                -4: "Warning: Two or three numbers separated by ':'.", \
                -3: "Warning: The negative numbers are too small.", \
                -2: "Warning: Above the total steps.", \
                -1: "Warning: Input error."})
    def _steps(self, steps):
        if isinstance(steps, Iterable):
            for step in [_-self._start+1 for _ in steps]:
                if steps > 0:
                    self.__steps.append(step)
                else:
                    self.__steps.append(step+self._total_steps+1)
        elif isinstance(steps, int):
            steps = steps - self._start + 1
            if steps > 0:
                self.__steps.append(steps)
            else:
                self.__steps.append(steps+self._total_steps+1)
    def _interval(self, steps):
        step_range = steps.split(':')
        if len(step_range) == 2:
            step_range[0] = int(step_range[0]) + 1
            step_range[1] = int(step_range[1]) + 2
        elif len(step_range) == 3:
            step_range[0] = int(step_range[0]) + 1
            step_range[1] = int(step_range[1])
            step_range[2] = int(step_range[2]) + 2
        if int(step_range[0]) < 0:
            step_range[0] += self._total_steps
            if len(step_range) == 2:
                step_range[1] += self._total_steps
            elif len(step_range) == 3:
                step_range[2] += self._total_steps
        else:
            step_range[0] -= self._start
            if len(step_range) == 2:
                step_range[1] -= self._start
            elif len(step_range) == 3:
                step_range[2] -= self._start
        steps = list(range(step_range[0], step_range[2], step_range[1])) if len(step_range) == 3 else list(range(step_range[0], step_range[1]))
        self.__steps.extend(steps)

class Step_LAMMPS(Step):
    def __init__(self, input_obj=None, put_flag=False, **kwargs):
        super().__init__(input_obj=input_obj, put_flag=put_flag, **kwargs)
    def _set_total_steps(self):
        from os import popen
        with popen(f"grep TIMESTEP {self._args.input} 2> /dev/null | wc") as command:
            self._total_steps = int(command.read().split()[0])

class Step_Gromacs(Step):
    def __init__(self, input_object=None, put_flag=False, **kwargs):
        super().__init__(input_object=input_object, put_flag=put_flag, **kwargs)
    def _set_total_steps(self):
        from os import popen
        with popen(f"grep step {self._args.input} | wc") as command:
            self._total_steps = int(command.read().split()[0])

class Atom(Fundamental):
    def __init__(self, input_obj=None, put_flag=False, **kwargs):
        super().__init__(input_obj=input_obj, **kwargs)
        del self._AtomStep
        self._atom_list = []
        self._index_dict = defaultdict(_none_factory)
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
            self._start = self._args.atom if hasattr(self._args, "atom") else 0
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
            for index, atom in enumerate(self._atom_list):
                self._index_dict[atom] = index
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
        self.args, self.start = software.args, software.args.atom if hasattr(software.args, "atom") else 1
        if hasattr(software, "atoms_info"):
            self.atoms_info = software.atoms_info
        if hasattr(software, "molecules"):
            self.molecules = software.molecules
        if hasattr(software, 'read_elements'):
            self.read_elements(software)
        elif software.elements != []:
            self.elements = software.elements
    @start.setter
    def start(self, start):
        if fullmatch(r"\-?\d+", str(start)) is not None:
            self._start = int(start)
        else:
            raise ValueError("The start number should be an integer.")
    @elements.setter
    def elements(self, elements):
        if isinstance(elements, list):
            self._elements = elements
        elif isinstance(elements, tuple):
            self._elements = list(elements)
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
            for atom, index in index_list.items():
                self._index_dict[atom] = index
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
        if fullmatch(r"\-?\d+", str(input_atom_list)) is not None:
            if int(input_atom_list) > self._atom_num:
                return -3
            self._atoms(input_atom_list)
        elif fullmatch(r"\-?\d+(_\-?\d+)*", str(input_atom_list)) is not None:
            self._atoms(input_atom_list.split('_'))
        elif isinstance(input_atom_list, str) and all([element in self._elements for element in input_atom_list.split("_")]):
            self._atom_kind(input_atom_list)
        elif ':' in input_atom_list:
            if fullmatch(r"\-?\d+:\-?\d+(:\-?\d+)?", str(input_atom_list)) is None:
                return -5
            sp = [int(atom)-self._start for atom in input_atom_list.split(':')]
            if (len(sp) == 3 and int(sp[0]) > int(sp[2])) or (len(sp) == 2 and int(sp[0]) > int(sp[1])):
                return -2
            if min(sp) < 0 and min(sp)+self._atom_num < 0:
                return -4
            elif max(sp) > self._atom_num:
                return -3
            else:
                self._interval(input_atom_list)
        elif "_except_" in str(input_atom_list):
            original_atoms_str, except_atoms_str = input_atom_list.split("_except_")
            original_atoms, except_atoms = Atom(self._args), Atom(self._args)
            original_atoms.elements, except_atoms.elements = self._elements, self._elements
            original_atoms.put(original_atoms_str); except_atoms.put(except_atoms_str)
            if len(original_atoms.get()) == 0 or len(except_atoms.get()) == 0:
                return -1
            elif len(set(original_atoms.get()) - set(except_atoms.get())) == 0:
                return -1
            for atom in list(set(original_atoms.get()) - set(except_atoms.get())):
                self._in_list_check(atom)
        elif fullmatch(r"[0-9A-Za-z_]+", str(input_atom_list)) is not None:
            if input_atom_list.lower() == "end":
                return 1
            elif input_atom_list.lower() == "all":
                self._atom_list = list(range(self._atom_num))
                return 2
            elif "mol" in str(input_atom_list).lower():
                if self._molecules.molecule_dictionary == {}:
                    if not hasattr(self._args, "molecules"):
                        return -6
                    elif self._args.molecules is not None and path.isfile(self._args.molecules):
                        self._molecules.read_molecules(self._args.molecules)
                        self.read_atoms_info(self._args.molecules)
                    else:
                        return -7
                return self._mol(input_atom_list)
            else:
                return -1
        elif isinstance(input_atom_list, Iterable):
            if not str(input_atom_list).strip():
                return -1
            elif fullmatch(r"\-?\d+", str(input_atom_list[0])) is not None:
                return self._atoms(input_atom_list)
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
        return defaultdict(int, {-8: "Warning: Can not match the molecule type.", \
                -7: "Warning: The molecule file does not exist.", \
                -6: "Warning: 'molecules' argument is not defined.", \
                -5: "Warning: Two or three numbers separated by ':'.", \
                -4: "Warning: The negative numbers are too small.", -3: "Warning: Above the total steps.", \
                -2: "Warning: The start number should be less than the end number.", -1: "Warning: Input error."})
    def _atoms(self, atoms):
        if isinstance(atoms, int) or (isinstance(atoms, str) and fullmatch(r"-?\d+(\.\d+)?", str(atoms)) is not None):
            atoms = int(atoms)
            if atoms >= 0:
                self._in_list_check(atoms-self._start)
            else:
                self._in_list_check(atoms+self._atom_num)
        elif isinstance(atoms, Iterable) and not isinstance(atoms, str):
            for atom in atoms:
                if int(atom) >= 0:
                    self._in_list_check(int(atom)-self._start)
                else:
                    self._in_list_check(int(atom)+self._atom_num-self._start)
        else:
            return -1
    def _mol(self, inp):
        inp = inp.replace("mol", "").replace("_", "", 1)
        if inp in [_.split('_')[0] for _ in self._molecules.molecule_kind]:
            for i, kind in enumerate(self._molecules.molecule_kind):
                if kind.split('_')[0] == inp:
                    for atom in self._molecules.molecule_dictionary[i]:
                        self._in_list_check(atom)
        elif len(inp.split('_')) > 1 and f"{inp.split('_')[0]}_{inp.split('_')[1]}" in self._molecules.molecule_kind:
            if len(inp.split('_')) == 2:
                for atom in self._molecules.molecule_dictionary[self._molecules.molecule_kind.index(inp)]:
                    self._in_list_check(atom)
            else:
                for element in inp.split("_")[2:]:
                    for atom in self._molecules.molecule_dictionary[self._molecules.molecule_kind.index(f"{inp.split('_')[0]}_{inp.split('_')[1]}")]:
                        if element == self._elements[atom]:
                            self._in_list_check(atom)
        elif len([element in self._elements for element in inp.split('_')[1:]]) > 0 and all([element in self._elements for element in inp.split('_')[1:]]):
            for mol, kind in enumerate(self._molecules.molecule_kind):
                if inp.split('_')[0] == kind.split('_')[0]:
                    for element in inp.split("_")[1:]:
                        for atom in self._molecules.molecule_dictionary[mol]:
                            if element == self._elements[atom]:
                                self._in_list_check(atom)
        else:
            return -8
        return 0
    def _interval(self, inp):
        atom_range = inp.split(':')
        if len(atom_range) == 2:
            atom_range[0] = int(atom_range[0])
            atom_range[1] = int(atom_range[1]) + 1
        elif len(atom_range) == 3:
            atom_range[0] = int(atom_range[0])
            atom_range[1] = int(atom_range[1])
            atom_range[2] = int(atom_range[2]) + 1
        if int(atom_range[0]) < 0:
            atom_range[0] += self._atom_num
            if len(atom_range) == 2:
                atom_range[1] += self._atom_num
            elif len(atom_range) == 3:
                atom_range[2] += self._atom_num
        else:
            atom_range[0] -= self._start
            if len(atom_range) == 2:
                atom_range[1] -= self._start
            elif len(atom_range) == 3:
                atom_range[2] -= self._start
        atom_list = list(range(atom_range[0], atom_range[1])) if len(atom_range) == 2 else list(range(atom_range[0], atom_range[2], atom_range[1]))
        self._atom_list.extend(atom_list)
    def _atom_kind(self, inp):
        atoms = [idx for idx, element in enumerate(self._elements) if element in inp.split('_')]
        if len(atoms) > 0:
            for atom in atoms:
                self._in_list_check(atom)
        else:
            return -1
    def _in_list_check(self, atom):
        if atom not in self._atom_list:
            self._atom_list.append(atom)
    def read_atoms_info(self, filepath=None):
        if filepath is None:
            if not hasattr(self._args, 'molecules'):
                raise FileNotFoundError("'molecules' file not found.")
            else:
                if path.isfile(self._args.molecules):
                    filepath = self._args.molecules
                else:
                    raise FileNotFoundError("'molecules' file not found.")
        with open(filepath) as read_file:
            tmp = []
            molecule = 0
            for row in read_file:
                row = row.replace('\n', '')
                sp = split(r"[^a-zA-Z0-9]", row)
                if "mol" in sp[0]:
                    tmp.append(sp[1])
                elif fullmatch(r"\d+(\,\d+)*", row) is not None:
                    for _ in [int(_) for _ in row.split(',')]:
                        self._atoms_info[_]["molecule"] = molecule
                    molecule += 1
    def read_elements(self, software, filepath=None):
        self._elements = software.read_elements() if software.elements == [] else software.elements
        if filepath is not None or (hasattr(self._args, 'elementfile') and self._args.elementfile is not None):
            if filepath is None:
                if path.isfile(self._args.elementfile):
                    filepath = self._args.elementfile
                else:
                    raise FileNotFoundError("'elementfile' file not found.")
            with open(filepath) as read_file:
                self._elements = list(self._elements)
                for idx, line in enumerate(read_file):
                    line, sp = line.replace('\n', ''), line.replace('\n', '').split()
                    if idx == 0 and fullmatch(r"[A-Za-z]+( [A-Za-z]+)*", line) is not None:
                        self._elements = sp
                        break
                    elif fullmatch(r"[A-Za-z0-9]+: [A-Za-z0-9_:,]+", line) is not None:
                        atoms = Atom(self._args)
                        atoms.elements = self._elements
                        for inp in sp[1].split(','):
                            err_return = atoms.put(inp)
                            if err_return:
                                print(f"line {idx+1}:", atoms.err_list[err_return], "Skip this line.")
                                break
                        else:
                            for atom in atoms.get():
                                self._elements[atom] = sp[0].replace(':', '')
                    elif not line.strip():
                        continue
                    else:
                        print(f"Warning: Input proper format in line {idx+1}. Skip this line.")
                self._elements = tuple(self._elements)
        self._atom_num = len(self._elements)

def grid(cen_pos: np.ndarray,
         mea_pos: np.ndarray,
         lattice: np.ndarray,
         wrap: np.ndarray):
    """
    cen_pos: (N_cen, 1, 1, 3)
    mea_pos: (1, N_mea, 1, 3)
    lattice: (3, 3)
    wrap:    (1, 1, N_wrap, 3)
    """

    # 直接使用 broadcasting 加總後乘 lattice
    cen_pos_new = np.matmul(cen_pos, lattice)  # (N_cen, 1, 1, 3)
    mea_pos_new = np.matmul(mea_pos + wrap, lattice)  # (1, N_mea, N_wrap, 3)

    return cen_pos_new, mea_pos_new

def distance_matrix_wrap(cen_pos: np.ndarray,
                         mea_pos: np.ndarray,
                         lattice: np.ndarray,
                         wrap: np.ndarray) -> np.ndarray:
    """
    cen_pos: (N_cen, 3)
    mea_pos: (N_mea, 3)
    lattice: (3, 3)
    wrap:    (N_wrap, 3)
    """

    # broadcasting reshape
    cen_pos_b = cen_pos[:, None, None, :]   # (N_cen, 1, 1, 3)
    mea_pos_b = mea_pos[None, :, None, :]   # (1, N_mea, 1, 3)
    wrap_b    = wrap[None, None, :, :]      # (1, 1, N_wrap, 3)

    # 套用 grid 計算座標
    cen, mea = grid(cen_pos_b, mea_pos_b, lattice, wrap_b)

    # 向量差
    vectors = cen - mea

    # 計算每一個向量的 L2 norm → (N_cen, N_mea, N_wrap)
    dist = np.linalg.norm(vectors, axis=3)

    return dist

class Molecule(Fundamental):
    def __init__(self, input_obj=None, **kwargs):
        super().__init__(input_obj=input_obj, **kwargs)
        del self._AtomStep
        self._mole_dict = defaultdict(list)
        self._mole_kind = []
        self._mole_pos_dict = defaultdict(list)
        self._bond_type = defaultdict(float)
        self._frac_flag, self._cart_flag = False, False
        self._atom_list = None
        if hasattr(input_obj, "bridge"):
            self.bridge = input_obj
        else:
            self.args = input_obj
        if hasattr(self._args, 'molecules') and self._args.molecules is not None and path.isfile(self._args.molecules):
            self.read_molecules(self._args.molecules)
    @property
    def bridge(self):
        pass
    @property
    def threshold(self):
        return self._bond_type
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
            self._atom_list = software.atoms.get()
            self._bond_type = software.molecules.threshold
    @threshold.setter
    def threshold(self, bond_type):
        if isinstance(bond_type, dict):
            if self._bond_type == {}:
                self._bond_type = bond_type
            else:
                for key, val in bond_type.items():
                    self._bond_type[key] = val
        else:
            raise Exception("'threshold' should be a dictionary.")
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
            if filepath is None:
                if hasattr(self._args, 'bond'):
                    if path.isfile(self._args.bond):
                        filepath = self._args.bond
                    else:
                        raise FileNotFoundError("'bond' file not found.")
                else:
                    raise ValueError("No 'bond' argument provided.")
            self.build_threshold()
            with open(filepath) as read_file:
                row = 0
                for line in read_file:
                    row += 1
                    line = line.replace('\n', '')
                    if fullmatch(r"[A-Za-z]{1,2}-[A-Za-z]{1,2} -?\d+\.?\d+", line) is not None:
                        sp = line.split()[0].split('-')
                        if sp[0] not in self._elements:
                            print(f"Warning: {sp[0]} in row {row} does not exist in this system. Skip this row.")
                        elif sp[1] not in self._elements:
                            print(f"Warning: {sp[1]} in row {row} does not exist in this system. Skip this row.")
                        elif float(line.split()[1]) < 0:
                            print(f"Warning: {line.split()[1]} in row {row} should be a positive number. Skip this row.")
                        else:
                            self._bond_type[f'{sp[0]}-{sp[1]}'] = self._bond_type[f'{sp[1]}-{sp[0]}'] = float(line.split()[1])
                    else:
                        print(f"Warning: Input proper format in row {row}. Skip this row.")
    def read_molecules(self, filepath=None):
        if self._mole_dict == {}:
            if not hasattr(self._args, 'molecules'):
                if filepath is None:
                    raise FileNotFoundError("'molecules' file not found.")
            elif filepath is None and self._args.molecules is not None:
                if path.isfile(self._args.molecules):
                    filepath = self._args.molecules
                else:
                    raise FileNotFoundError("'molecules' file not found.")
            with open(filepath) as read_file:
                mole_kind = []
                former_mole_kind = ""
                molecule = 0
                if self._atom_list is not None:
                    atom_dict = {atom: atom in self._atom_list for atom in range(len(self._elements))}
                for row in read_file:
                    row = row.replace('\n', '')
                    sp = split(r"[^a-zA-Z0-9]", row)
                    if "mol" in sp[0]:
                        mole_kind.append(sp[1])
                    elif fullmatch(r"\d+(\,\d+)*", row) is not None:
                        if mole_kind == []:
                            mole_kind.append("mol")
                        if mole_kind[-1] != former_mole_kind:
                            former_mole_kind = mole_kind[-1]
                            kind_counter = 1
                        else:
                            kind_counter += 1
                        if self._atom_list is not None:
                            for atom in sorted([int(_) for _ in row.split(',')]):
                                if atom_dict[atom-self._args.atom]:
                                    self._mole_dict[molecule].append(atom-self._args.atom)
                        else:
                            self._mole_dict[molecule] = sorted([int(atom)-self._args.atom for atom in row.split(',')])
                        self._mole_kind.append(mole_kind[-1]+f"_{kind_counter}")
                        molecule += 1
    def molecule_wrap(self):
        from kit.accelerate import image_shift, shift_to_origin
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
    def __init__(self, input_obj=None, atoms=None, steps=None, **kwargs):
        super().__init__(**kwargs)
        from kit.accelerate import c2f, f2c
        self._atoms = Atom(input_obj) if atoms is None else atoms
        self._steps = Step(input_obj) if steps is None else steps
        self._atom_dict = defaultdict(list)
        self._molecules = Molecule(input_obj)
        self._lattice = Lattice()
        self._cart_flag, self._frac_flag, self._fast_flag = False, False, False
        self.c2f, self.f2c = c2f, f2c
        if hasattr(input_obj, "bridge"):
            self.bridge = input_obj
        else:
            self.args = input_obj
        if hasattr(self._args, "molecules") and self._args.molecules is not None and path.isfile(self._args.molecules):
            self._molecules.read_molecules(self._args.molecules)
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
    def elements(self):
        return self._atoms.elements
    @property
    def atoms_info(self):
        return self._atoms.atoms_info
    @property
    def Lattice(self):
        return self._lattice
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
        self.Lattice = software.Lattice
    @atoms.setter
    def atoms(self, atoms):
        if isinstance(atoms, Atom):
            self._atoms = atoms
        else:
            raise ValueError(f"Only 'Atom' class can be imported (Imported type: {type(atoms)}).")
    @elements.setter
    def elements(self, elements):
        self._atoms.elements = elements
        self._molecules.elements = elements
    @atoms_info.setter
    def atoms_info(self, atoms_info):
        self._atoms.atoms_info = atoms_info
    @Lattice.setter
    def Lattice(self, lattice):
        if isinstance(lattice, Lattice):
            self._lattice = lattice
        else:
            raise ValueError(f"Only 'Lattice' class can be imported (Imported type: {type(lattice)}).")
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
    def fast_append(self, position, atom):
        self._tmp_pos.append(position)
    def put(self, position, atom):
        self._atom_dict[atom].append(position)
    @abstractmethod
    def split(self, atom_step, steps, atoms):
        if atom_step.__class__.__bases__[0] != AtomStep:
            raise ValueError("Only 'AtomStep_Trj' or 'AtomStep_Single_Point' can be imported.")
        atoms.sort()
        steps.sort()
        self.atoms = atoms
        self.steps = steps
        original_step, original_atom = atom_step.steps.get(), atom_step.atoms.get()
        original_step_dict = {step: True if step in original_step else False for step in steps.get()}
        original_atom_dict = {atom: True if atom in original_atom else False for atom in atoms.get()}
        positions = []
        step_dict, atom_dict = atom_step.steps.index_list, atom_step.atoms.index_list
        for step in steps.get():
            position = []
            if original_step_dict[step]:
                step_idx = step_dict[step]
                for atom in atoms.get():
                    if original_atom_dict[atom]:
                        atom_idx = atom_dict[atom]
                        position.append(atom_step.fractional_position[step_idx, atom_idx])
                positions.append(position)
        if atom_step.fractional_flag:
            self._frac_pos = array(positions)
            self._frac_flag = True
        else:
            self._cart_pos = array(positions)
            self._cart_flag = True
        if len(self._atoms.elements) != len(self._atoms.get()):
            self._atoms.elements = []
            for atom in self._atoms.get():
                self._atoms.elements.append(atom_step.atoms.elements[atom_dict[atom]])
            self._molecules.elements = self._atoms.elements
    @abstractmethod
    def lock(self):
        atom_list = self._atoms.get()
        if not hasattr(self._args, "molecules") and len(self._atoms.atoms_info) > 0 and "molecule" in self._atoms.atoms_info[atom_list[0]].keys():
            for key in atom_list:
                self._molecules.molecule_dictionary[self._atoms.atoms_info[key]["molecule"]].append(key)
                self._atoms.atoms_info[key].pop("molecule")
                self.build_molecule_position()
        if len(self._atoms.elements) and len(list(self._atoms.atoms_info.keys())) > 0 and hasattr(self._atoms.atoms_info[list(self._atoms.atoms_info.keys())[0]], "element"):
            for key in atom_list:
                self._atoms.elements.append(self._atoms.atoms_info[key]["element"])
                self._atoms.atoms_info[key].pop("element")
            self._molecules.elements = self._atoms.elements

        if self._frac_flag:
            self._frac_pos = []
            for atom in self._atoms.get():
                self._frac_pos.append(self._atom_dict[atom])
            self._frac_pos = array(self._frac_pos).transpose((1, 0, 2)) if array(self._frac_pos).ndim == 3 else array(self._frac_pos)
        if self._cart_flag:
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
        atom_dict = self._atoms.index_list
        for key, val in self._molecules.molecule_dictionary.items():
            if atom_dict[val[0]]:
                continue
            for atom in val:
                self._molecules.molecule_position[key].append(position[atom_dict[atom]])
            self._molecules.molecule_position[key] = array(self._molecules.molecule_position[key]).transpose((1, 0, 2)) if position.ndim == 3 else array(self._molecules.molecule_position[key])
        self._molecules.molecule_wrap()

class AtomStep_Trj(AtomStep):
    def __init__(self, input_obj=None, steps=None, atoms=None, **kwargs):
        super().__init__(input_obj=input_obj, atoms=atoms, steps=steps, **kwargs)
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
        assert self._cart_pos.shape[0] == len(steps), f"The first dimension of cartesian_position should be equal to the number of steps. ({self._cart_pos.shape[0]} != {len(steps)})"
        assert self._cart_pos.shape[1] == len(self._atoms.get()), f"The second dimension of cartesian_position should be equal to the number of atoms. ({self._cart_pos.shape[1]} != {len(self._atoms.get())})"
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
        assert self._frac_pos.shape[0] == len(steps), f"The first dimension of fractional_position should be equal to the number of steps. ({self._frac_pos.shape[0]} != {len(steps)})"
        assert self._frac_pos.shape[1] == len(self._atoms.get()), f"The second dimension of fractional_position should be equal to the number of atoms. ({self._frac_pos.shape[1]} != {len(self._atoms.get())})"
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
        assert self._fast_pos.shape[0] == len(steps), f"The first dimension of fast_position should be equal to the number of steps. ({self._fast_pos.shape[0]} != {len(steps)})"
        assert self._fast_pos.shape[1] == len(self._atoms.get()), f"The second dimension of fast_position should be equal to the number of atoms. ({self._fast_pos.shape[1]} != {len(self._atoms.get())})"
    @steps.setter
    def steps(self, steps):
        if isinstance(steps, Step):
            self._steps = steps
        else:
            raise ValueError("Only 'Step' class can be imported.")
    def split(self, atom_step, steps, atoms):
        super().split(atom_step, steps, atoms)
    def lock(self):
        super().lock()

class AtomStep_Single_Point(AtomStep):
    def __init__(self, input_obj=None, atoms=None, **kwargs):
        super().__init__(input_obj=input_obj, atoms=atoms, **kwargs)
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
        assert self._cart_pos.shape[0] == len(self._atoms.get()), f"The first dimension of cartesian_position should be equal to the number of atoms. ({self._cart_pos.shape[0]} != {len(self._atoms.get())})"
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
        assert self._frac_pos.shape[0] == len(self._atoms.get()), f"The first dimension of fractional_position should be equal to the number of atoms. ({self._frac_pos.shape[0]} != {len(self._atoms.get())})"
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
        assert self._fast_pos.shape[0] == len(self._atoms.get()), f"The first dimension of fast_position should be equal to the number of atoms. ({self._fast_pos.shape[0]} != {len(self._atoms.get())})"
    def split(self, atom_step, atoms):
        if atom_step.__class__.__bases__[0] != AtomStep:
            raise ValueError("Only 'AtomStep_Trj' or 'AtomStep_Single_Point' can be imported.")
        self.atoms = atoms
        original_atom = atom_step.atoms.get()
        if atom_step.fractional_flag:
            self._frac_pos = []
            for atom in atoms.get():
                if atom in original_atom:
                    atom_idx = atom_step.atoms.index_list[atom]
                    self._frac_pos.append(atom_step.fractional_position[atom_idx])
            self._frac_pos = array(self._frac_pos)
        else:
            self._cart_pos = []
            for atom in atoms.get():
                if atom in original_atom:
                    atom_idx = atom_step.atoms.index_list[atom]
                    self._cart_pos.append(atom_step.cartesian_position[atom_idx])
            self._cart_pos = array(self._cart_pos)
        if len(self._atoms.elements) != len(self._atoms.get()):
            self._atoms.elements = []
            for atom in self._atoms.get():
                self._atoms.elements.append(atom_step.atoms.elements[atom_step.atoms.index_list[atom]])
            self._molecules.elements = self._atoms.elements
    def lock(self):
        super().lock()
        if self._frac_flag:
            self._frac_pos = self._frac_pos[0]
        if self._cart_flag:
            self._cart_pos = self._cart_pos[0]
        if self._fast_flag:
            self._fast_pos = self._fast_pos[0]

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
    def __init__(self, input_obj=None, **kwargs):
        super().__init__(input_obj=input_obj, **kwargs)
        self.defaultdict = defaultdict
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
        from kit.accelerate import distance_matrix_cutoff, distance_matrix
        self.distance_matrix_cutoff, self.distance_matrix_func = distance_matrix_cutoff, distance_matrix
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
        return self._adj_list
    @property
    def molecules(self):
        return self._AtomStep.molecules
    @Fundamental.atom_step.setter
    def atom_step(self, atom_step):
        if isinstance(atom_step, AtomStep_Trj) or isinstance(atom_step, AtomStep_Single_Point):
            self._AtomStep = atom_step
            self.molecules = self._AtomStep.molecules
        else:
            raise ValueError("'atom_step' should be ' 'AtomStep_Trj' or 'AtomStep_Single_Point'.")
    @bridge.setter
    def bridge(self, software):
        self.args, self.elements, self.atom_step = software.args, software.elements, software.atom_step
        if hasattr(self._args, "molecules") and self._args.molecules is not None and self._AtomStep.molecules.molecule_dictionary != {}:
            if path.isfile(self._args.molecules):
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
        if isinstance(self._adj_list, list):
            self._adj_list = adjList
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
        atoms = self._AtomStep.atoms.get()
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
        atom_dict = self._AtomStep.atoms.index_list
        for mol_num, molecule in mole_dict.items():
            molecule.sort()
            if isinstance(mole_pos[mol_num], list) and mole_pos[mol_num] == []:
                if positions.ndim == 2:
                    for atom in molecule:
                        mole_pos[mol_num].append(positions[atom_dict[atom]])
                    mole_pos[mol_num] = array(mole_pos[mol_num])
                elif positions.ndim == 3:
                    for atom in molecule:
                        mole_pos[mol_num].append(positions[:, atom_dict[atom]])
                    mole_pos[mol_num] = array(mole_pos[mol_num]).transpose((1, 0, 2))
        if wrap:
            self._AtomStep.molecules.fractional_flag = True
            self._AtomStep.molecules.molecule_wrap()
        return mole_dict, mole_pos
    def build_graph(self, wrap=True, molecule=True):
        elements = self._AtomStep.atoms.elements
        self._adj_mat = np.zeros((len(elements), len(elements)))
        if self._AtomStep.molecules.threshold == {}:
            self._AtomStep.molecules.build_threshold()
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
                    self._dist_mat[idx] = self.distance_matrix_cutoff(position[np.newaxis, :], positions, self._AtomStep.lattice, self._period_images, max(self._AtomStep.molecules.threshold.values()), True)[0] if wrap else self.distance_matrix_func(position[np.newaxis, :], positions, self._AtomStep.lattice)
            self._adj_mat, self._adj_list = self.build_graph_from_distance_matrix(self._dist_mat)
    def build_distance_matrix(self, wrap=True):
        position = self._AtomStep.fractional_position if self._AtomStep.fractional_position.ndim == 2 else self._AtomStep.fractional_position[0]
        if wrap:
            return self.distance_matrix_cutoff(position, position, self._AtomStep.lattice, self._period_images, max(self._AtomStep.molecules.threshold.values()), True)[0]
        else:
            return self.distance_matrix_func(position, position, self._AtomStep.lattice)
    def molecule_distance_matrix(self):
        atoms_idx = self._AtomStep.atoms.index_list
        elements = self._AtomStep.atoms.elements
        for mole_num, mole_pos in self._AtomStep.molecules.molecule_position.items():
            if len(mole_pos) == 1:
                continue
            mole_dict = self._AtomStep.molecules.molecule_dictionary[mole_num]
            position = mole_pos if mole_pos.ndim == 2 else mole_pos[0]
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
        adj_list = defaultdict(list)
        atom_list = self._AtomStep.atoms.get()
        for i, j in zip(*np.where(adj_mat > 0)):
            atom_i, atom_j = atom_list[i], atom_list[j]
            adj_list[atom_i].append([atom_j, adj_mat[i, j]])
            adj_list[atom_j].append([atom_i, adj_mat[j, i]])
        return adj_list
    def list_to_matrix(self, adjacent_list=None):
        adj_list = adjacent_list if adjacent_list is not None and isinstance(adjacent_list, dict) else self._adj_list
        adj_mat = np.zeros((len(self._AtomStep.atoms.get()), len(self._AtomStep.atoms.get())))
        atoms_idx = self._AtomStep.atoms.index_list
        for cen_atom, mea_atoms in adj_list.items():
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
        atom_list = self._AtomStep.atoms.get()
        searched.append(sorted([atom_list[atom] for atom in mole_list]))
        if atom2 not in searched[0]:
            mole_list = self.walking(atom2_idx, graph)
            searched.append(sorted([atom_list[atom] for atom in mole_list]))
        if molecule is not None:
            tmp = [self._AtomStep.atoms.index_list[atom] for atom in molecule if atom not in searched[0] and atom not in searched[1]]
            while(tmp != []):
                atom = tmp.pop(0)
                mole_list = self.walking(atom, graph)
                searched.append(sorted([atom_list[atom] for atom in mole_list]))
                for atom in mole_list:
                    if atom in tmp:
                        tmp.remove(atom)
        return searched
    @staticmethod
    def walking(index, graph):
        new, mole_list = [index], []
        while(new != []):
            atom = new.pop()
            if atom not in mole_list:
                mole_list.append(atom)
            connected_atoms = list(np.where(graph[atom])[0])
            new.extend([connected_atom for connected_atom in connected_atoms if connected_atom not in mole_list])
        return mole_list
    @staticmethod
    def walking_list(index, graph):
        new, mole_list = [index], []
        while(new != []):
            atom = new.pop()
            if atom not in mole_list:
                mole_list.append(atom)
            connected_atoms = list(atom2 for atom2 in graph[atom].keys() if graph[atom][atom2])
            new.extend([connected_atom for connected_atom in connected_atoms if connected_atom not in mole_list])
        return mole_list
