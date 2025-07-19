from re import compile
from collections.abc import Iterable
from numpy import array, tile, zeros, linalg, pi, ndarray
import numpy as np
from kit.fundamental import Position, Periodic, AtomStep, AtomStep_Single_Point, AtomStep_Trj, Graph
from kit.accelerate import angle, dihedral_angle

class Wrap_Base(Periodic):
    def __init__(self, input_obj=None):
        self._directions = []
        self._period_images = array([[0, 0, 1], [1, 0, 1], [-1, 0, 1],
                                    [0, 1, 1], [0, -1, 1], [1, 1, 1],
                                    [1, -1, 1], [-1, 1, 1], [-1, -1, 1],
                                    [0, 0, 0], [1, 0, 0], [-1, 0, 0],
                                    [0, 1, 0], [0, -1, 0], [1, 1, 0],
                                    [1, -1, 0], [-1, 1, 0], [-1, -1, 0],
                                    [0, 0, -1], [1, 0, -1], [-1, 0, -1],
                                    [0, 1, -1], [0, -1, -1], [1, 1, -1],
                                    [1, -1, -1], [-1, 1, -1], [-1, -1, -1]])
        self._period_images_length = len(self._period_images)
        if hasattr(input_obj, "bridge"):
            self.args, self.atom_step = input_obj.args, input_obj.atom_step
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
    def center_atom_step(self):
        cen_pos = self._cen_pos if self._cen_pos.ndim == 2 else self._cen_pos.transpose((1, 0, 2))
        if isinstance(self._cen_atom_step, Iterable):
            for position, atom_step in zip(cen_pos, self._cen_atom_step):
                atom_step.fractional_position = position
        elif self._cen_atom_step.__class__.__bases__[0] == AtomStep:
            self._cen_atom_step.fractional_position = cen_pos
        return self._cen_atom_step
    @property
    def center_atom(self):
        if isinstance(self._cen_atom_step, Iterable):
            atoms = []
            for atom_step in self._cen_atom_step:
                atoms.append(atom_step.atoms)
        else:
            atoms = self._cen_atom_step.atoms
        return atoms
    @property
    def center_position(self):
        return self._cen_pos if self._cen_pos.ndim == 2 else self._cen_pos.transpose((1, 0, 2))
    @property
    def measure_atom_step(self):
        mea_pos = self._mea_pos if self._mea_pos.ndim == 2 else self._mea_pos.transpose((1, 0, 2))
        if isinstance(self._mea_atom_step, Iterable):
            for position, atom_step in zip(mea_pos, self._mea_atom_step):
                atom_step.fractional_position = position
        elif self._mea_atom_step.__class__.__bases__[0] == AtomStep:
            self._mea_atom_step.fractional_position = mea_pos
        return self._mea_atom_step
    @property
    def measure_atom(self):
        if isinstance(self._mea_atom_step, Iterable):
            atoms = []
            for atom_step in self._mea_atom_step:
                atoms.append(atom_step.atoms)
        else:
            atoms = self._mea_atom_step.atoms
        return atoms
    @property
    def measure_position(self):
        return self._mea_pos if self._mea_pos.ndim == 2 else self._mea_pos.transpose((1, 0, 2))
    @property
    def directions(self):
        return self._directions
    @atom_step.setter
    def atom_step(self, atom_step):
        if atom_step.__class__.__bases__[0] == AtomStep or isinstance(atom_step, Iterable):
            self._AtomStep = atom_step
        else:
            raise ValueError("Only 'AtomStep_Trj' or 'AtomStep_Single_Point' or 'Iterable' can be imported.")
    @bridge.setter
    def bridge(self, software):
        if isinstance(self._AtomStep, Iterable) and self._AtomStep[0].__class__.__bases__[0] == AtomStep:
            self.center_atom_step = software.atom_step[0]
            self.measure_atom_step = software.atom_step[1:]
        elif self._AtomStep.__class__.__bases__[0] == AtomStep:
            self._cen_atom_step = (AtomStep_Single_Point() if isinstance(self._AtomStep, AtomStep_Single_Point) else AtomStep_Trj())
            self._mea_atom_step = (AtomStep_Single_Point() if isinstance(self._AtomStep, AtomStep_Single_Point) else AtomStep_Trj())
            self._mea_atom_step.atoms.elements = self._cen_atom_step.atoms.elements = software.atom_step.atoms.elements
            atom_list = software.atom_step.atom_list
            self._cen_atom_step.atoms.put(atom_list[0])
            self._mea_atom_step.atoms.put(atom_list[1:])
            if isinstance(self._AtomStep, AtomStep_Trj):
                self._cen_atom_step.steps = software.steps
                self._mea_atom_step.steps = software.steps
    @period_images.setter
    def period_images(self, period_images):
        self._period_images = array(period_images)
        self._period_images_length = len(self._period_images)
    @directions.setter
    def directions(self, directions):
        self._directions = directions
    @center_atom_step.setter
    def center_atom_step(self, cen_atom_step):
        if isinstance(cen_atom_step, Iterable) and cen_atom_step[0].__class__.__bases__[0] == AtomStep:
            self._cen_atom_step = cen_atom_step
            self._cen_pos = []
            for cen_atom_step_iter in self._cen_atom_step:
                self._cen_pos.append(cen_atom_step_iter.fractional_position)
        elif cen_atom_step.__class__.__bases__[0] == AtomStep:
            self._cen_atom_step = cen_atom_step
            self._cen_pos = self._cen_atom_step.fractional_position
        else:
            raise ValueError("Only 'Iterable', 'AtomStep_Trj', or 'AtomStep_Single_Point' can be imported.")
    @center_position.setter
    def center_position(self, cen_pos):
        if not isinstance(cen_pos, ndarray):
            cen_pos = array(cen_pos)
        if isinstance(self._cen_atom_step, Iterable):
            if len(self._cen_atom_step) != len(cen_pos):
                raise ValueError(f"The number of positions ({len(cen_pos)} positions) should be equal to the number of atom_step ({len(self._cen_atom_step)} atom_step).")
            else:
                for atom_step, position in zip(self._cen_atom_step, cen_pos):
                    atom_step.fractional_position = position
                self._cen_pos = cen_pos
        else:
            self._cen_atom_step.fractional_position = cen_pos
            self._cen_pos = cen_pos[:, np.newaxis, :] if cen_pos.ndim == 2 else cen_pos.transpose((1, 0, 2))
    @measure_atom_step.setter
    def measure_atom_step(self, meaAtomStep):
        if isinstance(meaAtomStep, Iterable) and meaAtomStep[0].__class__.__bases__[0] == AtomStep:
            self._mea_atom_step = meaAtomStep
            self._mea_pos = []
            for MeaAtomStep in self._mea_atom_step:
                self._mea_pos.append(MeaAtomStep.fractional_position)
        elif meaAtomStep.__class__.__bases__[0] == AtomStep:
            self._mea_atom_step = meaAtomStep
            self._mea_pos = self._mea_atom_step.fractional_position
        else:
            raise ValueError("Only 'Iterable', 'AtomStep_Trj', or 'AtomStep_Single_Point' can be imported.")
    @measure_position.setter
    def measure_position(self, mea_pos):
        if not isinstance(mea_pos, ndarray):
            mea_pos = array(mea_pos)
        if isinstance(self._mea_atom_step, Iterable):
            if len(self._mea_atom_step) != len(mea_pos):
                raise ValueError(f"The number of positions ({len(mea_pos)} positions) should be equal to the number of atom_step ({len(self._mea_atom_step)} atom_step).")
            else:
                for atom_step, position in zip(self._mea_atom_step, mea_pos):
                    atom_step.fractional_position = position
                self._mea_pos = mea_pos if mea_pos.ndim == 2 else mea_pos.transpose((1, 0, 2))
        else:
            self._mea_atom_step.fractional_position = mea_pos
            self._mea_pos = mea_pos[:, np.newaxis, :] if mea_pos.ndim == 2 else mea_pos.transpose((1, 0, 2))
    def adjust(self):
        if self._cen_pos.ndim == 2:
            self._cen_pos[0] += self._directions[0]
            self._mea_pos[:][0] += self._directions[1:]
        elif self._cen_pos.ndim == 1:
            self._cen_pos += self._directions[0]
            self._mea_pos[:] += self._directions[1:]

class Wrap(Wrap_Base):
    def __init__(self, input_obj=None):
        from kit.accelerate import rev_image_shift, image_shift
        Wrap_Base.__init__(self, input_obj)
        self.rev_image_shift = rev_image_shift
        self.image_shift = image_shift
    def manual(self):
        self.adjust()
        self._cen_pos = self.image_shift(self._cen_pos, self._AtomStep.lattice)
        for idx, mea_pos in enumerate(self._mea_pos):
            self._mea_pos[idx] = self.image_shift(mea_pos, self._AtomStep.lattice)
    def auto(self, fractional_distance=0.75):
        for idx, mea_pos in enumerate(self._mea_pos):
            direction = ((mea_pos[0] - self._cen_pos[0][0]) > fractional_distance)
            mea_pos[0][direction] -= 1
            direction = ((mea_pos[0] - self._cen_pos[0][0]) < -fractional_distance)
            mea_pos[0][direction] += 1
            self._mea_pos[idx] = self.image_shift(mea_pos, self._AtomStep.lattice)
        self._cen_pos = self.image_shift(self._cen_pos, self._AtomStep.lattice)
    def rev_auto(self, fractional_distance=0.75):
        for idx, mea_pos in enumerate(self._mea_pos):
            direction = ((mea_pos[-1] - self._cen_pos[-1][0]) > fractional_distance)
            mea_pos[-1][direction] -= 1
            direction = ((mea_pos[-1] - self._cen_pos[-1][0]) < -fractional_distance)
            mea_pos[-1][direction] += 1
            self._mea_pos[idx] = self.rev_image_shift(mea_pos, self._AtomStep.lattice)
        self._cen_pos = self.rev_image_shift(self._cen_pos, self._AtomStep.lattice)
    def sprinkler(self):
        def expand(position, wrap_length):
            if position.ndim == 2:
                return tile(position, wrap_length).reshape(position.shape[0], wrap_length, 3)
            elif position.ndim == 3:
                return tile(position, wrap_length).reshape(position.shape[0], position.shape[1], wrap_length, 3)
        if isinstance(self._cen_pos, list):
            positions = []
            for position in self._cen_pos:
                positions.append(expand(position, self._period_images_length))
            self._cen_pos = positions
        elif isinstance(self._cen_pos, ndarray):
            self._cen_pos = expand(self._cen_pos, self._period_images_length)
        if isinstance(self._mea_pos, list):
            positions = []
            for position in self._mea_pos:
                positions.append(expand(position, self._period_images_length))
            self._mea_pos = positions
        elif isinstance(self._mea_pos, ndarray):
            self._mea_pos = expand(self._mea_pos, self._period_images_length)
    def benchmark(self, benchmark, fractional_distance=0.75):
        for idx, mea_pos in enumerate(self._mea_pos):
            direction = ((mea_pos[benchmark] - self._cen_pos[benchmark][0]) > fractional_distance)
            mea_pos[benchmark][direction] -= 1
            direction = ((mea_pos[benchmark] - self._cen_pos[benchmark][0]) < -fractional_distance)
            mea_pos[benchmark][direction] += 1
            self._mea_pos[idx][benchmark:] = self.image_shift(mea_pos[benchmark:], self._AtomStep.lattice)
            if benchmark != 0:
                self._mea_pos[idx][:(benchmark+1)] = self.rev_image_shift(mea_pos[:(benchmark+1)], self._AtomStep.lattice)

class Distance(Wrap):
    def __init__(self, input_obj=None):
        Wrap.__init__(self, input_obj)
        self._directions = [[0, 0, 0], [0, 0, 0]]
    def distance(self, center=None, measure=None, axis=1):
        if center is None:
            center = self._cen_pos
        if measure is None:
            measure = self._mea_pos
        if center.ndim == 2 and measure.ndim == 2:
            return linalg.norm(measure-center, axis=axis)
        elif center.ndim == 3 and measure.ndim == 3:
            distances = []
            for center_pos, measure_pos in zip(center, measure):
                distances.append(linalg.norm(measure_pos-center_pos, axis=axis))
            return array(distances)
        elif center.ndim != measure.ndim:
            raise ValueError(f"The dimension of central and measure position should be equal (center: {center.ndim} dim, measure: {measure.ndim} dim).")
        elif center.ndim > 3 or center.ndim < 2:
            raise ValueError("The dimension should be 2 or 3 for central position.")
        elif measure.ndim > 3 or measure.ndim < 2:
            raise ValueError("The dimension should be 3 or 4 for measure position.")

class Mobility(Periodic):
    def __init__(self):
        Periodic.__init__(self)
    @property
    def bridge(self):
        pass
    @bridge.setter
    def bridge(self, software):
        self.args, self.atom_step = software.args, software.atom_step
        self._lattice = self.atom_step.lattice
    def benchmark(self, benchmark=0):
        self.__distance = Distance()
        self.__distance.atom_step = self._AtomStep
        fraPos = self._AtomStep.fractional_position
        fraShape = self._AtomStep.fractional_position.shape
        carPos = zeros(fraShape)
        carPos[:, benchmark:] = self.__distance.image_shift(fraPos[:, benchmark:])
        carPos[:, :(benchmark+1)] = self.__distance.image_shift(fraPos[:, :(benchmark+1)])
        self._AtomStep.cartesian_position = carPos
        benchmark_pos = tile(fraPos[benchmark], (fraShape[1], 1)).reshape((fraShape[1], fraShape[0], 27))
        self.__benchmark_pos = self._AtomStep.f2c(benchmark_pos, self._lattice)
    def mobility(self):
        distances = []
        for idx, carPos in enumerate(self._AtomStep.cartesian_position):
            distances.append(self.__distance.distance(carPos, self.__benchmark_pos[idx]))
        return array(distances).T

class Angle(Wrap):
    def __init__(self, input_obj=None):
        Wrap.__init__(self, input_obj)
        self._directions = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    def angle(self, center=None, measure_1=None, measure_2=None):
        if center is None:
            center = self._cen_pos
        if measure_1 is None:
            measure_1 = self._mea_pos[0]
        if measure_2 is None:
            measure_2 = self._mea_pos[1]
        degreeFlag = (False if self._args.radian else True)
        if center.ndim == 2:
            return angle(measure_1-center, measure_2-center, degreeFlag)
        elif center.ndim == 3:
            angles = []
            for center_pos, measure_1_pos, measure_2_pos in zip(center, measure_1, measure_2):
                angles.append(angle(measure_1_pos-center_pos, measure_2_pos-center_pos, degreeFlag))
            return array(angles)
        elif center.ndim != measure_1.ndim or center.ndim != measure_2.ndim:
            raise ValueError(f"The dimension of central and measure position should be equal (center: {center.ndim} dim, side_1: {measure_1.ndim} dim, side_2: {measure_2.ndim} dim).")
        elif center.ndim > 3 or center.ndim < 2:
            raise ValueError("The dimension should be 2 or 3 for central position.")
        elif measure_1.ndim > 3 or measure_1.ndim < 2:
            raise ValueError("The dimension should be 3 or 4 for measure position.")

class Dihedral_Angle(Wrap):
    def __init__(self, input_obj=None):
        Wrap.__init__(self, input_obj)
        self._directions = [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]]
    def dihedral_angle(self, center=None, measure_1=None, measure_2=None, measure_3=None):
        if center is None:
            center = self._cen_pos
        if measure_1 is None:
            measure_1 = self._mea_pos[0]
        if measure_2 is None:
            measure_2 = self._mea_pos[1]
        if measure_3 is None:
            measure_3 = self._mea_pos[2]
        degreeFlag = (False if self._args.radian else True)
        if center.ndim == 2:
            return dihedral_angle(measure_1-center, measure_2-center, measure_3-center, degreeFlag)
        elif center.ndim == 3:
            dihedral_angles = []
            for center_pos, measure_1_pos, measure_2_pos, measure_3Pos in zip(center, measure_1, measure_2, measure_3):
                dihedral_angles.append(dihedral_angle(measure_1_pos-center_pos, measure_2_pos-center_pos, measure_3Pos-center_pos, degreeFlag))
            return array(dihedral_angles)
        elif center.ndim != measure_1.ndim or center.ndim != measure_2.ndim or center.ndim != measure_3.ndim:
            raise ValueError(f"The dimension of central and measure position should be equal (center: {center.ndim} dim, collinear: {measure_1.ndim} dim, side_1: {measure_2.ndim} dim, side_2: {measure_3.ndim} dim).")
        elif center.ndim > 3 or center.ndim < 2:
            raise ValueError("The dimension should be 2 or 3 for central position.")
        elif measure_1.ndim > 3 or measure_1.ndim < 2:
            raise ValueError("The dimension should be 3 or 4 for measure position.")

class Smooth:
    @property
    def factor(self):
        return self.__factor
    @factor.setter
    def factor(self, factor):
        if factor.lower() == "vh":
            self.__factor = 0.997
        elif factor.lower() == "h":
            self.__factor = 0.99
        elif factor.lower() == "m":
            self.__factor = 0.8
        elif factor.lower() == "s":
            self.__factor = 0.6
        elif compile(r"^\d\.\d+$").match(str(factor)) is None:
            print("Warning: Input a number between 0 and 1.")
            self.__factor = -1
        elif float(factor) > 1:
            print("Warning: This number should not be larger than 1.")
            self.__factor = -1
        else:
            self.__factor = float(factor)
    def smooth(self, curve, factor):
        from kit.accelerate import EMA
        unsmoothed_curves = array(curve).T if array(curve).shape[0] > array(curve).shape[1] * 10 else array(curve)
        self.factor = factor
        return EMA(unsmoothed_curves, self.__factor)

class Coordination:
    def __init__(self):
        self.__MoleWeight, self.__LiquidDensity, self.__density = None, None, None
    @property
    def radius(self):
        return self.__radius
    @property
    def g(self):
        return self.__g
    @property
    def molecule_weight(self):
        return self.__MoleWeight
    @property
    def liquid_density(self):
        return self.__LiquidDensity
    @property
    def density(self):
        return self.__density
    @radius.setter
    def radius(self, radius):
        self.__radius = array(radius)
    @g.setter
    def g(self, g):
        self.__g = array(g)
    @density.setter
    def density(self, density):
        self.__density = density
    @molecule_weight.setter
    def molecule_weight(self, MoleWeight):
        self.__MoleWeight = MoleWeight
    @liquid_density.setter
    def liquid_density(self, LiquidDensity):
        self.__LiquidDensity = LiquidDensity
    def cumulative_function(self):
        if (self.__MoleWeight is not None) and (self.__LiquidDensity is not None):
            self.__density = self.__LiquidDensity/self.__MoleWeight
        elif self.__density is None:
            raise Exception("Provide information about molecule weight or liquid density.")
        from scipy.interpolate import interp1d
        from scipy.integrate import quad
        import scipy.constants as C
        def inner_function(radius, rho):
            integ = radius[0]
            InterpFunc = interp1d(radius, rho, kind="cubic")
            AccumuFunc = zeros(radius.shape)
            AccumuFunc[0] = integ
            for delta in range(1, len(radius)):
                tmp = quad(InterpFunc, radius[delta-1], radius[delta])[0]
                if tmp > 0:
                    integ += tmp
                AccumuFunc[delta] = integ
            return AccumuFunc
        if self.__g.ndim == 2:
            if self.__g.shape[0] > self.__g.shape[1]:
                self.__g = self.__g.T
            self.__rho = 4*pi*self.__density*10**(-24)*C.N_A*(self.__radius)**2*self.__g
            AccumuFunc = zeros(self.__rho.shape)
            for i in range(self.__g.shape[0]):
                AccumuFunc[i] = inner_function(self.__radius, self.__rho[i])
            return AccumuFunc.T
        elif self.__g.ndim == 1:
            self._rho = 4*pi*self.__density*10**(-24)*C.N_A*(self.__radius)**2*self.__g
            return inner_function(self.__radius, self.__rho)
    def coordination_number(self, coordNum):
        if (self.__MoleWeight is not None) and (self.__LiquidDensity is not None):
            self.__density = self.__LiquidDensity/self.__MoleWeight
        elif self.__density is None:
            raise Exception("Provide information about molecule weight or liquid density.")
        from scipy.interpolate import interp1d
        from scipy.integrate import quad
        import scipy.constants as C
        def inner_function(radius, rho, coordNum):
            integ = radius[0]
            InterpFunc = interp1d(radius, rho, kind="cubic")
            for delta in range(1, len(radius)):
                tmp = quad(InterpFunc, radius[delta-1], radius[delta])[0]
                if integ < coordNum and (integ+tmp) >= coordNum:
                    diff = coordNum - integ
                    step = np.linspace(radius[delta-1], radius[delta], 201)[1] - np.linspace(radius[delta-1], radius[delta], 201)[0]
                    for i in np.linspace(radius[delta-1], radius[delta], 201):
                        if quad(InterpFunc, radius[delta-1], i)[0] < diff and quad(InterpFunc, radius[delta-1], i+step)[0] >= diff:
                            return i+step
                if tmp > 0:
                    integ += tmp
        if self.__g.ndim == 2:
            coordNums = []
            if self.__g.shape[0] > self.__g.shape[1]:
                self.__g = self.__g.T
            self.__rho = 4*np.pi*self.__density*10**(-24)*C.N_A*(self.__radius)**2*self.__g
            for i in range(self.__g.shape[0]):
                coordNums.append(inner_function(self.__radius, self.__rho[i], coordNum))
            return coordNums
        elif self.__g.ndim == 1:
            self.__rho = 4*np.pi*self.__density*10**(-24)*C.N_A*(self.__radius**2)*self.__g
            return inner_function(self.__radius, self.__rho[i], coordNum)

class Solvent(Graph):
    def __init__(self, input_obj=None):
        Graph.__init__(self, input_obj)
        from collections import defaultdict
        self.defaultdict = defaultdict
        if hasattr(input_obj, "bridge"):
            self._lattice = input_obj.atom_step.lattice
            if hasattr(input_obj.args, "scaling"):
                self.search_scaling = input_obj.args.scaling
            else:
                self._search_scaling = 2.5
        else:
            self._search_scaling = 2.5
    @property
    def bridge(self):
        pass
    @property
    def search_scaling(self):
        return self._search_scaling
    @bridge.setter
    def bridge(self, software):
        Graph.bridge.fset(self, software)
        self._lattice = software.atom_step.lattice
        if hasattr(software.args, "scaling"):
            self.search_scaling = software.args.scaling
    @search_scaling.setter
    def search_scaling(self, value):
        if compile(r"^\d+\.?\d*$").match(str(value)) is not None:
            self._search_scaling = float(value)
        else:
            raise ValueError("The search scaling must be a positiove number.")
    def average_position(self, molecule_position, axis=1):
        avg_pos = {}
        for mol, position in molecule_position.items():
            if position is not None:
                avg_pos[mol] = np.average(position, axis=axis)
        return avg_pos

class Ligand(Solvent):
    def __init__(self, input_obj=None):
        Solvent.__init__(self, input_obj)
    @property
    def lattice_bond(self):
        return self._lattice_bond
    @property
    def cluster_graph(self):
        return self._cluster_adj_list
    @property
    def ligand_molecules(self):
        return self._ligand_molecules
    @property
    def cluster_atoms(self):
        return self._cluster_atoms
    @property
    def cluster_position(self):
        return self._cluster_pos
    @lattice_bond.setter
    def lattice_bond(self, lattice_bond):
        if compile(r"^\d+\.?\d*$").match(str(lattice_bond)) is not None:
            self._lattice_bond = lattice_bond
        else:
            raise ValueError("The lattice bond should be a positive number.")
    def ligand_search(self, cation=None, anion=None):
        axis = 0 if self._AtomStep.molecules.molecule_position[0].ndim == 2 else 1
        avg_pos = self.average_position(self._AtomStep.molecules.molecule_position, axis=axis)
        
        self._cen_mole, self._cen_avg_pos, self._mea_mole, self._mea_avg_pos = [], [], [], []
        cen_atom, mea_atom = self._cen_atom_step.atoms.get(), self._mea_atom_step.atoms.get()
        read_atoms = []
        for atom in range(len(self._AtomStep.elements)):
            if atom in read_atoms:
                continue
            for key, val in self._AtomStep.molecules.molecule_dictionary.items():
                if atom in cen_atom and atom in val and key not in self._cen_mole:
                    self._cen_avg_pos.append(avg_pos[key])
                    self._cen_mole.append(key)
                    read_atoms.extend(val)
                if atom in mea_atom and atom in val and key not in self._mea_mole:
                    self._mea_avg_pos.append(avg_pos[key])
                    self._mea_mole.append(key)
                    read_atoms.extend(val)
        
        if self._args.mode.lower() == "g":
            cation_atoms, anion_atoms = cation.get(), anion.get()
            self._cation_moles, self._anion_moles, self._solvent_moles = [], [], []
            for mole_num, mole_dict in self._AtomStep.molecules.molecule_dictionary.items():
                for atom in mole_dict:
                    if atom in cation_atoms:
                        self._cation_moles.append(mole_num)
                    elif atom in anion_atoms:   
                        self._anion_moles.append(mole_num)
                    else:
                        self._solvent_moles.append(mole_num)
                    break
        
        self._cen_avg_pos, self._mea_avg_pos = array(self._cen_avg_pos), array(self._mea_avg_pos)
        
        if self._AtomStep.fractional_position.ndim == 2:
            self._build_cluster_graph(self._AtomStep.molecules.molecule_position, self._cen_avg_pos, self._mea_avg_pos)
            self._ligand_molecules, self._cluster_atoms, self._cluster_pos = self._graph_to_cluster_position(self._cluster_adj_list, self._AtomStep.molecules.molecule_dictionary, self._AtomStep.molecules.molecule_position)
        elif self._AtomStep.fractional_position.ndim == 3:
            self._ligand_molecules, self._cluster_atoms, self._cluster_pos = self.defaultdict(list), self.defaultdict(list), self.defaultdict(list)
            for t, (cen_avg_pos, mea_avg_pos) in enumerate(zip(self._cen_avg_pos.transpose((1, 0, 2)), self._mea_avg_pos.transpose((1, 0, 2)))):
                mole_pos = {mol: position[t] for mol, position in self._AtomStep.molecules.molecule_position.items()}
                self._build_cluster_graph(mole_pos, cen_avg_pos, mea_avg_pos)
                ligand, cluster_atoms, cluster_pos = self._graph_to_cluster_position(self._cluster_adj_list, self._AtomStep.molecules.molecule_dictionary, mole_pos)
                for mol in ligand.keys():
                    self._ligand_molecules[mol].append(ligand[mol]); self._cluster_atoms[mol].append(cluster_atoms[mol]); self._cluster_pos[mol].append(cluster_pos[mol])
    def _build_cluster_graph(self, molecule_position, center_avg_molecules, measure_avg_molecules):
        self._cluster_adj_list = self.defaultdict(dict)
        threshold = self._lattice_bond * self._search_scaling
        lattice = self._AtomStep.lattice
        dist_mat = self.distance_matrix_wrap(center_avg_molecules, measure_avg_molecules, lattice, self._period_images)
        cluster_candidates = self.defaultdict(list)
        for i, j, _ in zip(*np.where(dist_mat < threshold)):
            cluster_candidates[self._cen_mole[i]].append(self._mea_mole[j])
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
    def build_neighbors(self, cluster_graph):
        ligand_molecules = self.defaultdict(list)
        for i in cluster_graph.keys():
            for j in cluster_graph[i].keys():
                for _ in range(cluster_graph[i][j]):
                    ligand_molecules[i].append(j)
        return ligand_molecules
    def _graph_to_cluster_position(self, cluster_graph, molecule_dictionary, molecule_position):
        ligand_molecules, cluster_pos, cluster_atom  = self.defaultdict(list), self.defaultdict(list), self.defaultdict(list)
        recorded = []
        
        for i in self._cen_mole:
            ligand_molecules[i].append(i)
            for j in cluster_graph[i].keys():
                ligand_molecules[i].append(j)
                cen_pos, mea_pos = molecule_position[i], molecule_position[j]
                dist_mat = self.distance_matrix_wrap(cen_pos, mea_pos, self._AtomStep.lattice, self._period_images)
                if i not in recorded:
                    cluster_atom[i].extend(molecule_dictionary[i])
                    cluster_pos[i].extend(cen_pos)
                    recorded.append(i)
                cluster_atom[i].extend(molecule_dictionary[j])
                for _, jj, kk in zip(*np.where(dist_mat < self._lattice_bond)):
                    mea_pos[jj] += self._period_images[kk]
                    mea_pos = self.image_shift(mea_pos, self._AtomStep.lattice, benchmark=jj, cartesian=0)
                    break
                cluster_pos[i].extend(mea_pos)

        for mol, positions in cluster_pos.items():
            cluster_pos[mol] = self._AtomStep.f2c(array(positions), self.lattice) - self._AtomStep.f2c(array(positions[0]), self.lattice)
        
        return ligand_molecules, cluster_atom, cluster_pos
    def cluster_types(self):
        if not hasattr(self, "_ligand_type"):
            if self._args.mode.lower() == "a":
                self._search_ligand_type(self._atom_type)
            elif self._args.mode.lower() == "m" or self._args.mode.lower() == "g":
                self._search_ligand_type(self._molecule_type)
        return (self._mole_kind_sorted, self._ligand_type)
    def _atom_type(self, ligand_molecules):
        ligand_type = self.defaultdict(list)
        for mol, ligands in ligand_molecules.items():
            ligand_type[mol] = [0]*len(self._mole_kind_sorted)
            for ligand in ligands:
                ligand_type[mol][self._mole_kind_sorted.index(self._AtomStep.molecules.molecule_kind[ligand].split("_")[0])] += 1
        return ligand_type
    def _molecule_type(self, ligand_molecules):
        ligand_type = self.defaultdict(list)
        for mol, ligands in ligand_molecules.items():
            ligand_type[mol] = [0]*len(self._mole_kind_sorted)
            for ligand in list(set(ligands)):
                ligand_type[mol][self._mole_kind_sorted.index(self._AtomStep.molecules.molecule_kind[ligand].split("_")[0])] += 1
        return ligand_type
    def _search_ligand_type(self, func=None):
        mole_kind = [kind.split('_')[0] for kind in self._AtomStep.molecules.molecule_kind]
        self._mole_kind_sorted = sorted(list(set(mole_kind)))
        if self._AtomStep.fractional_position.ndim == 2:
            self._ligand_type = func(self._ligand_molecules)
        elif self._AtomStep.fractional_position.ndim == 3:
            self._ligand_type = self.defaultdict(list)
            for t in range(self._AtomStep.fractional_position.shape[0]):
                ligand_molecules = {mol: ligands[t] for mol, ligands in self._ligand_molecules.items()}
                ligand_type = func(ligand_molecules)
                for mol, type in ligand_type.items():
                    self._ligand_type[mol].append(type)

class RDF(Wrap_Base):
    def __init__(self):
        Wrap_Base.__init__(self)
    @property
    def interval(self):
        return self._interval
    @property
    def rdf(self):
        return self._rdf
    @property
    def coordination_number(self):
        return self._coordNum
    def run(self, nbins=500, r_min=0.0, r_max=10.0):
        from kit.accelerate import distance_matrix_wrap, shift_to_origin
        lattice = self._AtomStep.lattice
        self._interval = np.linspace(r_min, r_max, nbins + 1)
        shell_volumes = 4.0 / 3.0 * np.pi * (self._interval[1:]**3 - self._interval[:-1]**3)
        if lattice.ndim == 2:
            volume = np.cross(lattice[0], lattice[1]).dot(lattice[2])
        else:
            volume = zeros(lattice.shape)
            for idx, lattice in enumerate(lattice):
                volume[idx] = np.cross(lattice[0], lattice[1]).dot(lattice[2])
        self._rdf, self._coordNum = zeros(nbins), zeros(nbins)
        self._cen_pos = shift_to_origin(self._cen_pos)
        self._mea_pos = shift_to_origin(self._mea_pos)
        cen_atom, mea_atom = self._cen_atom_step.atoms.get(), self._mea_atom_step.atoms.get()
        for fs, (cen_pos, mea_pos) in enumerate(zip(self._cen_pos, self._mea_pos)):
            if lattice.ndim == 2:
                rho = float(len(cen_atom)) / volume
                dij = distance_matrix_wrap(cen_pos, mea_pos, lattice, self._period_images)
            elif lattice.ndim == 3:
                rho = float(len(cen_atom)) / volume[fs]
                dij = distance_matrix_wrap(cen_pos, mea_pos, lattice[fs], self._period_images)
            hist = np.histogram(dij, bins=nbins, range=(r_min, r_max), density=False)[0]
            self._rdf += hist * 1.0 / rho
            self._coordNum += np.cumsum(hist)
        self._rdf /=  shell_volumes * len(self._cen_pos) * len(mea_atom)
        self._coordNum /= len(self._cen_pos) * len(mea_atom)
