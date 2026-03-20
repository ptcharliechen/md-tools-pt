from re import fullmatch
from copy import deepcopy
from collections.abc import Iterable
from numpy import array, tile, zeros, linalg, pi, ndarray
import numpy as np
from kit.fundamental import Position, Periodic, AtomStep, AtomStep_Single_Point, AtomStep_Trj, Graph
from kit.accelerate import angle, dihedral_angle

class Wrap_Base(Periodic):
    def __init__(self, input_obj=None, **kwargs):
        super().__init__(input_obj=input_obj, **kwargs)
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
        self.args, self.atom_step = software.args, software.atom_step
        if isinstance(self._AtomStep, Iterable) and self._AtomStep[0].__class__.__bases__[0] == AtomStep:
            self.center_atom_step = software.atom_step[0]
            self.measure_atom_step = software.atom_step[1:]
        elif self._AtomStep.__class__.__bases__[0] == AtomStep:
            from kit.fundamental import Atom
            self._cen_atom_step = (AtomStep_Single_Point() if isinstance(self._AtomStep, AtomStep_Single_Point) else AtomStep_Trj())
            self._mea_atom_step = (AtomStep_Single_Point() if isinstance(self._AtomStep, AtomStep_Single_Point) else AtomStep_Trj())
            self._mea_atom_step.atoms.elements = self._cen_atom_step.atoms.elements = software.atom_step.atoms.elements
            atom_list = software.atom_step.atom_list
            cen_atom, mea_atom = Atom(software), Atom(software)
            cen_atom.put(atom_list[0]); mea_atom.put(atom_list[1:])
            self._cen_atom_step.atoms = cen_atom
            self._mea_atom_step.atoms = mea_atom
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
            self.center_position = deepcopy(self._cen_atom_step.fractional_position)
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
                self._cen_pos = array(cen_pos)[np.newaxis, :, :] if cen_pos.ndim == 2 else cen_pos.transpose((1, 0, 2))
        else:
            self._cen_atom_step.fractional_position = cen_pos
            self._cen_pos = cen_pos[np.newaxis, :, :] if cen_pos.ndim == 2 else cen_pos.transpose((1, 0, 2))
    @measure_atom_step.setter
    def measure_atom_step(self, mea_atom_step):
        if isinstance(mea_atom_step, Iterable) and mea_atom_step[0].__class__.__bases__[0] == AtomStep:
            self._mea_atom_step = mea_atom_step
            self._mea_pos = []
            for MeaAtomStep in self._mea_atom_step:
                self._mea_pos.append(MeaAtomStep.fractional_position)
        elif mea_atom_step.__class__.__bases__[0] == AtomStep:
            self._mea_atom_step = mea_atom_step
            self.measure_position = deepcopy(self._mea_atom_step.fractional_position)
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
                self._mea_pos = array(mea_pos)[np.newaxis, :, :] if mea_pos.ndim == 2 else mea_pos.transpose((1, 0, 2))
        else:
            self._mea_atom_step.fractional_position = mea_pos
            self._mea_pos = mea_pos[np.newaxis, :, :] if mea_pos.ndim == 2 else mea_pos.transpose((1, 0, 2))
    def adjust(self):
        if self._cen_pos.ndim == 3:
            self._cen_pos[:] += self._directions[0]
            self._mea_pos[:][0] += self._directions[1:]
        elif self._cen_pos.ndim == 2:
            self._cen_pos[0] += self._directions[0]
            self._mea_pos[:][0] += self._directions[1:]
        elif self._cen_pos.ndim == 1:
            self._cen_pos += self._directions[0]
            self._mea_pos[:] += self._directions[1:]

class Wrap(Wrap_Base):
    def __init__(self, input_obj=None, **kwargs):
        from kit.accelerate import rev_image_shift, image_shift
        super().__init__(input_obj=input_obj, **kwargs)
        self.rev_image_shift = rev_image_shift
        self.image_shift = image_shift
    def manual(self):
        self.adjust()
        for idx, mea_pos in enumerate(self._mea_pos):
            self._mea_pos[idx] = self.image_shift(mea_pos, self._AtomStep.lattice)
        for idx, cen_pos in enumerate(self._cen_pos):
            self._cen_pos[idx] = self.image_shift(cen_pos, self._AtomStep.lattice)
    def auto(self, fractional_distance=0.75):
        for idx, mea_pos in enumerate(self._mea_pos):
            direction = ((mea_pos[0] - self._cen_pos[0][0]) > fractional_distance)
            mea_pos[0][direction] -= 1
            direction = ((mea_pos[0] - self._cen_pos[0][0]) < -fractional_distance)
            mea_pos[0][direction] += 1
            self._mea_pos[idx, :] = self.image_shift(mea_pos, self._AtomStep.lattice)
        for idx, cen_pos in enumerate(self._cen_pos):
            self._cen_pos[idx] = self.image_shift(cen_pos, self._AtomStep.lattice)
    def rev_auto(self, fractional_distance=0.75):
        for idx, mea_pos in enumerate(self._mea_pos):
            direction = ((mea_pos[-1] - self._cen_pos[-1][0]) > fractional_distance)
            mea_pos[-1][direction] -= 1
            direction = ((mea_pos[-1] - self._cen_pos[-1][0]) < -fractional_distance)
            mea_pos[-1][direction] += 1
            self._mea_pos[idx, :] = self.rev_image_shift(mea_pos, self._AtomStep.lattice)
        for idx, cen_pos in enumerate(self._cen_pos):
            self._cen_pos[idx] = self.rev_image_shift(cen_pos, self._AtomStep.lattice)
        #self._cen_pos = self.rev_image_shift(self._cen_pos, self._AtomStep.lattice)
    def sprinkler(self):
        def expand(position, wrap_length):
            if position.ndim == 2:
                return tile(position, wrap_length).reshape(position.shape[0], wrap_length, 3)
            elif position.ndim == 3:
                return tile(position, wrap_length).reshape(position.shape[0], position.shape[1], wrap_length, 3)
        for idx, cen_pos in enumerate(self._cen_pos):
            self._cen_pos[idx, :] = self.image_shift(cen_pos, self._AtomStep.lattice, cartesian=0)
        for idx, mea_pos in enumerate(self._mea_pos):
            self._mea_pos[idx, :] = self.image_shift(mea_pos, self._AtomStep.lattice, cartesian=0)
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
            self._mea_pos = positions + self._period_images
        elif isinstance(self._mea_pos, ndarray):
            self._mea_pos = expand(self._mea_pos, self._period_images_length) + self._period_images
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
        super().__init__(input_obj=input_obj)
        self._directions = [[0, 0, 0], [0, 0, 0]]
    def distance(self, center=None, measure=None, axis=1):
        if center is None:
            center = self._cen_pos
        if measure is None:
            measure = self._mea_pos
        if center.ndim != measure.ndim:
            raise ValueError(f"The dimension of central and measure position should be equal (center: {center.ndim} dim, measure: {measure.ndim} dim).")
        elif center.ndim == 2 and measure.ndim == 2:
            return linalg.norm(measure-center, axis=axis)
        elif center.ndim == 3 and measure.ndim == 3:
            distances = []
            for center_pos, measure_pos in zip(center, measure):
                distances.append(linalg.norm(measure_pos-center_pos, axis=axis))
            return array(distances)
        elif center.ndim == 4 and measure.ndim == 4:
            from kit.accelerate import f2c_acc
            distances = []
            for center_atom_pos in center:
                distance_arrays = []
                for measure_atom_pos in measure:
                    distance_array = []
                    for center_pos, measure_pos in zip(f2c_acc(center_atom_pos, self._AtomStep.lattice), f2c_acc(measure_atom_pos, self._AtomStep.lattice)):
                        distance_array.append(linalg.norm(measure_pos-center_pos, axis=axis))
                    distance_arrays.append(distance_array)
                distances.append(distance_arrays)
            return array(distances)
        elif center.ndim > 4 or center.ndim < 2:
            raise ValueError("The dimension should be 2, 3, or 4 for central position.")
        elif measure.ndim > 4 or measure.ndim < 2:
            raise ValueError("The dimension should be 2, 3, or 4 for measure position.")

class Mobility(Periodic):
    def __init__(self, input_obj=None, **kwargs):
        super().__init__(input_obj=input_obj, **kwargs)
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
            center = self._cen_pos[0]
        if measure_1 is None:
            measure_1 = self._mea_pos[0]
        if measure_2 is None:
            measure_2 = self._mea_pos[1]
        degree_flag = (False if self._args.radian else True)
        if center.ndim == 2:
            return angle(measure_1-center, measure_2-center, degree_flag)
        elif center.ndim == 3:
            angles = []
            for center_pos, measure_1_pos, measure_2_pos in zip(center, measure_1, measure_2):
                angles.append(angle(measure_1_pos-center_pos, measure_2_pos-center_pos, degree_flag))
            return array(angles)
        elif center.ndim != measure_1.ndim or center.ndim != measure_2.ndim:
            raise ValueError(f"The dimension of central and measure position should be equal (center: {center.ndim} dim, side_1: {measure_1.ndim} dim, side_2: {measure_2.ndim} dim).")
        elif center.ndim > 3 or center.ndim < 2:
            raise ValueError("The dimension should be 2 or 3 for central position.")
        elif measure_1.ndim > 3 or measure_1.ndim < 2:
            raise ValueError("The dimension should be 3 or 4 for measure position.")

class Dihedral_Angle(Wrap):
    def __init__(self, input_obj=None, **kwargs):
        super().__init__(input_obj=input_obj, **kwargs)
        self._directions = [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]]
    def dihedral_angle(self, center=None, measure_1=None, measure_2=None, measure_3=None):
        if center is None:
            center = self._cen_pos[0]
        if measure_1 is None:
            measure_1 = self._mea_pos[0]
        if measure_2 is None:
            measure_2 = self._mea_pos[1]
        if measure_3 is None:
            measure_3 = self._mea_pos[2]
        degree_flag = (False if self._args.radian else True)
        if center.ndim == 2:
            return dihedral_angle(measure_1-center, measure_2-center, measure_3-center, degree_flag)
        elif center.ndim == 3:
            dihedral_angles = []
            for center_pos, measure_1_pos, measure_2_pos, measure_3Pos in zip(center, measure_1, measure_2, measure_3):
                dihedral_angles.append(dihedral_angle(measure_1_pos-center_pos, measure_2_pos-center_pos, measure_3Pos-center_pos, degree_flag))
            return array(dihedral_angles)
        elif center.ndim != measure_1.ndim or center.ndim != measure_2.ndim or center.ndim != measure_3.ndim:
            raise ValueError(f"The dimension of central and measure position should be equal (center: {center.ndim} dim, collinear: {measure_1.ndim} dim, side_1: {measure_2.ndim} dim, side_2: {measure_3.ndim} dim).")
        elif center.ndim > 3 or center.ndim < 2:
            raise ValueError("The dimension should be 2 or 3 for central position.")
        elif measure_1.ndim > 3 or measure_1.ndim < 2:
            raise ValueError("The dimension should be 3 or 4 for measure position.")

class Direction(Position):
    def __init__(self, input_obj=None, **kwargs):
        super().__init__(input_obj=input_obj, **kwargs)
    @property
    def bridge(self):
        pass
    @bridge.setter
    def bridge(self, software):
        self.args, self.atom_step = software.args, software.atom_step
    @property
    def x_direction(self):
        self._wrap()
        return self._AtomStep.cartesian_position[:, :, 0]
    @property
    def y_direction(self):
        self._wrap()
        return self._AtomStep.cartesian_position[:, :, 1]
    @property
    def z_direction(self):
        self._wrap()
        return self._AtomStep.cartesian_position[:, :, 2]
    def _wrap(self):
        from kit.accelerate import image_shift
        self._AtomStep.cartesian_position = image_shift(self._AtomStep.fractional_position, self._AtomStep.lattice)

class Smooth:
    @property
    def factor(self):
        return self.__factor
    @property
    def curve(self):
        return self.__curves
    @factor.setter
    def factor(self, factor):
        if factor.lower() == "vh":
            self.__factor = 0.997
        elif factor.lower() == "h":
            self.__factor = 0.99
        elif factor.lower() == "m":
            self.__factor = 0.8
        elif factor.lower() == "l":
            self.__factor = 0.6
        elif fullmatch(r"\d\.\d+", str(factor)) is None:
            print("Warning: Input a number between 0 and 1.")
            self.__factor = -1
        elif float(factor) > 1:
            print("Warning: This number should not be larger than 1.")
            self.__factor = -1
        else:
            self.__factor = float(factor)
    @curve.setter
    def curve(self, curve):
        curve = array(curve)
        self.__unsmoothed_curves = curve if curve.ndim == 1 else curve.T
    def smooth(self):
        from kit.accelerate import EMA
        self.__curves = EMA(self.__unsmoothed_curves, self.__factor)

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
    def __init__(self, input_obj=None, **kwargs):
        super().__init__(input_obj=input_obj, **kwargs)
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
        if fullmatch(r"\d+\.?\d*", str(value)) is not None:
            self._search_scaling = float(value)
        else:
            raise ValueError("The search scaling must be a positiove number.")
    def average_position(self, molecule_position, axis=1):
        avg_pos = {}
        for mol, position in molecule_position.items():
            if position is not None:
                avg_pos[mol] = np.average(position, axis=axis)
        return avg_pos

class Reaction(Solvent):
    def __init__(self, input_obj=None, **kwargs):
        super().__init__(input_obj=input_obj, **kwargs)
    @property
    def bond_info(self):
        return self._bond_info
    def bond_breaking(self, t):
        mole_dict = self._fragment_mole_dict
        atom_list = self._AtomStep.atoms.get()
        for mol, position in self._fragment_mole_pos.items():
            dist_mat = self.distance_matrix_func(position[t], position[t], self.lattice)
            for idx_i, cen_atom in enumerate(mole_dict[mol]):
                cen_atom_idx = self._AtomStep.atoms.index_list[cen_atom]
                for mea_atom_idx in np.where(self._adj_mat[cen_atom_idx] > 0)[0]:
                    mea_atom = atom_list[mea_atom_idx]
                    idx_j = mole_dict[mol].index(atom_list[mea_atom_idx])
                    if idx_i < idx_j:
                        continue
                    if dist_mat[idx_i, idx_j] > self._adj_mat[cen_atom_idx, mea_atom_idx]:
                        self._adj_mat = self.delete_bond(cen_atom, mea_atom, self._adj_mat)
                        atom_min, atom_max = min(cen_atom, mea_atom), max(cen_atom, mea_atom)
                        fragments = self.connect(cen_atom, mea_atom, self._adj_mat)
                        position_min = self._AtomStep.f2c(position[t][idx_i], self.lattice) if cen_atom == atom_min else self._AtomStep.f2c(position[t][idx_j], self.lattice)
                        position_max = self._AtomStep.f2c(position[t][idx_i], self.lattice) if cen_atom == atom_max else self._AtomStep.f2c(position[t][idx_j], self.lattice)
                        bond_type = 1 if len(fragments) == 1 else 2
                        atom_min_idx, atom_max_idx = self._AtomStep.atoms.index_list[atom_min], self._AtomStep.atoms.index_list[atom_max]
                        if bond_type == 2:
                            elements_min = "".join([self._AtomStep.atoms.elements[idx] for idx in self.walking(atom_min_idx, self._adj_mat)])
                            elements_max = "".join([self._AtomStep.atoms.elements[idx] for idx in self.walking(atom_max_idx, self._adj_mat)])
                        distance = np.linalg.norm(position_min - position_max)
                        info = [t-2+self._args.step, atom_min, atom_max, position_min, position_max, distance]
                        if bond_type == 2:
                            info.append(elements_min)
                            info.append(elements_max)
                        if not any([bond_info[0] == info[0] and bond_info[1] == info[1] and bond_info[2] == info[2] for bond_info in self._bond_info[bond_type]]):
                            self._bond_info[bond_type].append(info)
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
        dist_mat, period_images = self.distance_matrix_cutoff(avg_pos_arr, avg_pos_arr, self.lattice, self._period_images, max(self._AtomStep.molecules.threshold.values())*self._search_scaling, True, 8)
        avg_pos_keys = list(avg_pos.keys())
        elements = self._AtomStep.atoms.elements
        for mole_idx_i, mole_idx_j in zip(*np.where((dist_mat < max(self._AtomStep.molecules.threshold.values())*self._search_scaling) & (dist_mat > 0.01))):
            mole_i, mole_j = avg_pos_keys[mole_idx_i], avg_pos_keys[mole_idx_j]
            frac_pos_i, frac_pos_j = mole_pos_dict[mole_i][t], mole_pos_dict[mole_j][t] + self._period_images[period_images[mole_idx_i, mole_idx_j]]
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
                    atom_idx = self._AtomStep.atoms.index_list[atom_min]
                    info = [t-2+self._args.step, atom_min, atom_max, position_min, position_max, distance, "".join([self._AtomStep.atoms.elements[idx] for idx in self.walking(atom_idx, self._adj_mat)])]
                    if not any([bond_info[0] == info[0] and bond_info[1] == info[1] and bond_info[2] == info[2] for bond_info in self._bond_info[3]]):
                        self._bond_info[3].append(info)
                        if not np.all(self._period_images[period_images[mole_idx_i, mole_idx_j]] == np.array([0, 0, 0])):
                            recorded.add((mole_j, period_images[mole_idx_i, mole_idx_j], idx_j))
        for info in list(recorded):
            mole_pos_dict[info[0]][t][info[2]] += self._period_images[info[1]]
            mole_pos_dict[info[0]][t] = self.image_shift(mole_pos_dict[info[0]][t], self.lattice, benchmark=info[2], cartesian=0)
            mole_pos_dict[info[0]][t:] = self.image_shift(mole_pos_dict[info[0]][t:], self.lattice, cartesian=0)
        for mol_idx in sorted(list(mole_dict.keys())):
            mol, position = mole_dict[mol_idx], mole_pos_dict[mol_idx]
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
                    atom_idx = self._AtomStep.atoms.index_list[atom_min]
                    info = [t-2+self._args.step, atom_min, atom_max, position_min, position_max, distance, ""]
                    if not any([bond_info[0] == info[0] and bond_info[1] == info[1] and bond_info[2] == info[2] for bond_info in self._bond_info[3]]):
                        self._bond_info[3].append(info)
        self.update_graph()
    def run(self):
        self._bond_info = {1: [], 2: [], 3: []}
        for t in range(len(self._AtomStep.steps.get())):
            self.bond_breaking(t)
            self.bond_creating(t)
        self._bond_info[1].sort(key=lambda x: x[0]); self._bond_info[2].sort(key=lambda x: x[0]); self._bond_info[3].sort(key=lambda x: x[0])

class Ligand(Solvent, Wrap_Base):
    def __init__(self, input_obj=None, **kwargs):
        super().__init__(input_obj=input_obj, **kwargs)
        from kit.accelerate import distance_matrix_wrap, distance_matrix_cutoff
        self.distance_matrix_wrap, self.distance_matrix_cutoff = distance_matrix_wrap, distance_matrix_cutoff
    @property
    def dative_bond(self):
        return self._dative_bond
    @dative_bond.setter
    def dative_bond(self, dative_bond):
        if fullmatch(r"\d+\.?\d*", str(dative_bond)) is not None:
            self._dative_bond = float(dative_bond)
        else:
            raise ValueError("The dative bond must be a positiove number.")
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
    def run(self, cation=None, anion=None):
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
        
        if self._args.mode.lower() == "c":
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
        threshold = self._dative_bond * self._search_scaling
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
                s = np.sum(dist_mat < self._dative_bond)
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
                if cluster_graph[i][j] > 1:
                    for _ in range(cluster_graph[i][j] - 1):
                        ligand_molecules[i].append(j)
                cen_pos, mea_pos = molecule_position[i], molecule_position[j]
                dist_mat = self.distance_matrix_wrap(cen_pos, mea_pos, self._AtomStep.lattice, self._period_images)
                if i not in recorded:
                    cluster_atom[i].extend(molecule_dictionary[i])
                    cluster_pos[i].extend(cen_pos)
                    recorded.append(i)
                cluster_atom[i].extend(molecule_dictionary[j])
                for _, jj, kk in zip(*np.where(dist_mat < self._dative_bond)):
                    mea_pos[jj] += self._period_images[kk]
                    mea_pos = self.image_shift(mea_pos, self._AtomStep.lattice, benchmark=jj, cartesian=0)
                    break
                cluster_pos[i].extend(mea_pos)

        for mol, positions in cluster_pos.items():
            cluster_pos[mol] = self._AtomStep.f2c(array(positions), self.lattice) - self._AtomStep.f2c(array(positions[0]), self.lattice)
        return ligand_molecules, cluster_atom, cluster_pos
    def cluster_types(self):
        mole_kind = [kind.split('_')[0] for kind in self._AtomStep.molecules.molecule_kind]
        self._mole_kind_sorted = sorted(list(set(mole_kind)))
        if not hasattr(self, "_ligand_type"):
            if self._args.mode.lower() == "a":
                self._search_ligand_type(self._atom_type)
            elif self._args.mode.lower() == "m" or self._args.mode.lower() == "c":
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
        if self._AtomStep.fractional_position.ndim == 2:
            self._ligand_type = func(self._ligand_molecules)
        elif self._AtomStep.fractional_position.ndim == 3:
            self._ligand_type = self.defaultdict(list)
            for t in range(self._AtomStep.fractional_position.shape[0]):
                ligand_molecules = {mol: ligands[t] for mol, ligands in self._ligand_molecules.items()}
                ligand_type = func(ligand_molecules)
                for mol, type in ligand_type.items():
                    self._ligand_type[mol].append(type)

class HCE(Ligand):
    def __init__(self, input_obj=None, **kwargs):
        super().__init__(input_obj=input_obj, **kwargs)
    def _build_cluster_graph(self, molecule_position, center_avg_molecules, measure_avg_molecules):
        self._cluster_adj_list = self.defaultdict(dict)
        threshold = self._dative_bond * self._search_scaling
        lattice = self._AtomStep.lattice
        cluster_candidates = self.defaultdict(list)
        for cen_mol_idx, cen_avg_mole in zip(self._cen_mole, center_avg_molecules):
            if cen_mol_idx in self._cation_moles:
                mea_mole = self._mea_mole
                mea_pos = measure_avg_molecules
            else:
                mea_mole, mea_pos = [], []
                for mea_mol_idx, pos in zip(self._mea_mole, measure_avg_molecules):
                    if mea_mol_idx in self._cation_moles:
                        mea_mole.append(mea_mol_idx)
                        mea_pos.append(pos)
            
            dist_mat = self.distance_matrix_wrap(cen_avg_mole[np.newaxis, :], array(mea_pos), lattice, self._period_images)
            for _, j, _ in zip(*np.where(dist_mat < threshold)):
                cluster_candidates[cen_mol_idx].append(mea_mole[j])
        
        for cen_mol_idx, mea_moles in cluster_candidates.items():
            cen_pos = molecule_position[cen_mol_idx]
            for mea_mol_id in mea_moles:
                if mea_mol_id == cen_mol_idx:
                    continue
                mea_pos = molecule_position[mea_mol_id]
                dist_mat = self.distance_matrix_wrap(cen_pos, mea_pos, lattice, self._period_images)
                s = np.sum(dist_mat < self._dative_bond)
                if s:
                    self._cluster_adj_list[mea_mol_id][cen_mol_idx] = self._cluster_adj_list[cen_mol_idx][mea_mol_id] = s
    def _graph_to_cluster_position(self, cluster_graph, molecule_dictionary, molecule_position):
        cluster_pos, cluster_atoms, ligand_molecules = self.defaultdict(list), self.defaultdict(list), self.defaultdict(list)
        
        recorded_cations = []
        for anion_mole in self._anion_moles:
            cluster_moles = self.walking_list(anion_mole, cluster_graph)
            cluster_avg_pos = []
            for cluster_mole in cluster_moles:
                cluster_atoms[anion_mole].extend(molecule_dictionary[cluster_mole])
                ligand_molecules[anion_mole].append(cluster_mole)
                if cluster_mole in cluster_graph[anion_mole].keys() and cluster_graph[anion_mole][cluster_mole] > 1:
                    for _ in range(cluster_graph[anion_mole][cluster_mole] - 1):
                        ligand_molecules[anion_mole].append(cluster_mole)
                mea_pos, mea_avg_pos = molecule_position[cluster_mole], np.average(molecule_position[cluster_mole], axis=0)
                if cluster_avg_pos != []:
                    dist_mat = self.distance_matrix_wrap(array(cluster_avg_pos), mea_avg_pos[np.newaxis, :], self._AtomStep.lattice, self._period_images)
                    _, _, k = np.unravel_index(np.argmin(dist_mat), dist_mat.shape)
                    mea_pos[0] += self._period_images[k]
                    mea_pos = self.image_shift(mea_pos, self._AtomStep.lattice, cartesian=0)
                cluster_pos[anion_mole].extend(list(mea_pos))
                cluster_avg_pos.append(mea_avg_pos)
                if cluster_mole != anion_mole and cluster_mole in cluster_atoms.keys():
                    del cluster_atoms[cluster_mole], cluster_pos[cluster_mole]
                if cluster_mole in self._cation_moles:
                    recorded_cations.append(cluster_mole)
                    ligand_molecules[cluster_mole] = ligand_molecules[anion_mole]
            for mol in ligand_molecules[anion_mole]:
                if mol in self._cation_moles:
                    ligand_molecules[mol] = ligand_molecules[anion_mole]
        
        for cation_mole in self._cation_moles:
            if cation_mole in recorded_cations:
                continue
            cluster_moles = self.walking_list(cation_mole, cluster_graph)
            cen_pos = molecule_position[cation_mole]
            for cluster_mole in cluster_moles:
                ligand_molecules[cation_mole].append(cluster_mole)
                cluster_atoms[cation_mole].extend(molecule_dictionary[cluster_mole])
                if cluster_mole in cluster_graph[cation_mole].keys():
                    for _ in range(cluster_graph[cation_mole][cluster_mole]):
                        ligand_molecules[cation_mole].append(cluster_mole)
                mea_pos = molecule_position[cluster_mole]
                dist_mat = self.distance_matrix_wrap(cen_pos, mea_pos, self._AtomStep.lattice, self._period_images)
                for _, jj, kk in zip(*np.where(dist_mat < self._dative_bond)):
                    mea_pos[jj] += self._period_images[kk]
                    mea_pos = self.image_shift(mea_pos, self._AtomStep.lattice, benchmark=jj, cartesian=0)
                    break
                cluster_pos[cation_mole].extend(list(mea_pos))
        for mol, positions in cluster_pos.items():
            cluster_pos[mol] = list(self._AtomStep.f2c(array(positions), self.lattice) - self._AtomStep.f2c(np.average(positions, axis=0), self.lattice))
        return ligand_molecules, cluster_atoms, cluster_pos
    def walking_list(self, index, graph):
        new, mole_list = [index], []
        while(new != []):
            atom = new.pop(0)
            if atom not in mole_list:
                mole_list.append(atom)
            if atom in self._solvent_moles:
                continue
            elif atom in self._anion_moles:
                connected_atom = list(atom2 for atom2 in graph[atom].keys() if graph[atom][atom2] and atom2 in self._cation_moles)
            else:
                connected_atom = list(atom2 for atom2 in graph[atom].keys() if graph[atom][atom2])
            new.extend([atom for atom in connected_atom if atom not in mole_list])
        return mole_list

class RDF(Wrap_Base):
    def __init__(self, input_obj=None, **kwargs):
        super().__init__(input_obj=input_obj, **kwargs)
    @property
    def interval(self):
        return self._interval
    @property
    def rdf(self):
        return self._rdf
    @property
    def coordination_number(self):
        return self._coord_num
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
        self._rdf, self._coord_num = zeros(nbins), zeros(nbins)
        self._cen_pos = shift_to_origin(array(self._cen_pos))
        self._mea_pos = shift_to_origin(array(self._mea_pos))
        cen_pos = self._cen_pos.transpose(1, 0, 2)[0]
        cen_atom, mea_atom = self._cen_atom_step.atoms.get(), self._mea_atom_step.atoms.get()
        for fs, (cen_pos, mea_pos) in enumerate(zip(self._cen_pos.transpose(1, 0, 2), self._mea_pos.transpose(1, 0, 2))):
            if lattice.ndim == 2:
                rho = float(len(cen_atom)) / volume
                dij = distance_matrix_wrap(cen_pos, mea_pos, lattice, self._period_images)
            elif lattice.ndim == 3:
                rho = float(len(cen_atom)) / volume[fs]
                dij = distance_matrix_wrap(cen_pos, mea_pos, lattice[fs], self._period_images)
            hist = np.histogram(dij, bins=nbins, range=(r_min, r_max), density=False)[0]
            self._rdf += hist * 1.0 / rho
            self._coord_num += np.cumsum(hist)
        self._rdf /=  shell_volumes * self._cen_pos.shape[1] * len(mea_atom)
        self._coord_num /= self._cen_pos.shape[1] * len(mea_atom)
        return self._rdf, self._coord_num

class Bond(Wrap_Base):
    def __init__(self, input_obj=None, **kwargs):
        super().__init__(input_obj=input_obj, **kwargs)
    def run(self):
        from kit.accelerate import shift_to_origin, distance_matrix_cutoff
        self._AtomStep.fractional_position = shift_to_origin(self._AtomStep.fractional_position)
        avg_bond, std_arr = 0, []
        cen_atom_list = self._cen_atom_step.atoms.get()
        for time, (cen_pos, mea_pos) in enumerate(zip(self._cen_atom_step.fractional_position, self._mea_atom_step.fractional_position)):
            dist_mat = distance_matrix_cutoff(cen_pos, mea_pos, self._AtomStep.lattice[time] if self._AtomStep.lattice.ndim == 3 else self._AtomStep.lattice, self._period_images, max(self._AtomStep.molecules.threshold.values()), True)[0]
            counter = {key: 0 for key in self._cen_atom_step.atoms.get()}
            for i, j in zip(*np.where((dist_mat < max(self._AtomStep.molecules.threshold.values())) & (dist_mat > 0.01))):
                element_i = self._cen_atom_step.atoms.elements[i]
                element_j = self._mea_atom_step.atoms.elements[j]
                bond_type = f"{element_i}-{element_j}"
                if dist_mat[i, j] < self._AtomStep.molecules.threshold[bond_type]:
                    avg_bond += 1
                    counter[cen_atom_list[i]] += 1
            std_arr.extend(list(counter.values()))
        avg_bond /= self._cen_pos.shape[0] * self._cen_pos.shape[1]
        return avg_bond, np.std(std_arr)

class Bond_graph(Graph, Wrap_Base):
    def __init__(self, input_obj=None, **kwargs):
        super().__init__(input_obj=input_obj, **kwargs)
    def run(self):
        from kit.accelerate import shift_to_origin
        self._AtomStep.fractional_position = shift_to_origin(self._AtomStep.fractional_position)
        cluster = self.defaultdict(list)
        atom_list = self._AtomStep.atoms.get()
        for time, position in enumerate(self._AtomStep.fractional_position):
            dist_mat = self.distance_matrix_cutoff(position, position, self._AtomStep.lattice[time] if self._AtomStep.lattice.ndim == 3 else self._AtomStep.lattice, self._period_images, max(self._AtomStep.molecules.threshold.values()), True)[0]
            for i, j in zip(*np.where(self._adj_mat > 0)):
                if i > j:
                    continue
                element_i = self._AtomStep.atoms.elements[i]
                element_j = self._AtomStep.atoms.elements[j]
                bond_type = f"{element_i}-{element_j}"
                if dist_mat[i, j] > self._AtomStep.molecules.threshold[bond_type]:
                    self._adj_mat[i, j] = self._adj_mat[j, i] = 0
            for i, j in zip(*np.where((dist_mat < max(self._AtomStep.molecules.threshold.values())) & (dist_mat > 0.01))):
                if i > j:
                    continue
                element_i = self._AtomStep.atoms.elements[i]
                element_j = self._AtomStep.atoms.elements[j]
                bond_type = f"{element_i}-{element_j}"
                if dist_mat[i, j] < self._AtomStep.molecules.threshold[bond_type] and self._adj_mat[i, j] == 0:
                    self._adj_mat[i, j] = self._adj_mat[j, i] = self._AtomStep.molecules.threshold[bond_type]
            recorded = []
            cen_atom_list = self._cen_atom_step.atoms.get()
            for cen_atom in cen_atom_list:
                if cen_atom in recorded:
                    continue
                connected_atom_list_idx = self.walking(self._AtomStep.atoms.index_list[cen_atom], self._adj_mat)
                connected_atom_list = [atom_list[atom] for atom in connected_atom_list_idx]
                cluster[cen_atom].append(connected_atom_list)
                recorded.append(cen_atom)
                for atom in connected_atom_list:
                    if atom in cen_atom_list and atom not in recorded:
                        cluster[atom].append(connected_atom_list)
                        recorded.append(atom)
        return cluster
    def walking(self, index, graph):
        new_idx, atom_list_idx = [index], []
        cen_atoms_idx = [self._AtomStep.atoms.index_list[atom] for atom in self._cen_atom_step.atoms.get()]
        mea_atoms_idx = [self._AtomStep.atoms.index_list[atom] for atom in self._mea_atom_step.atoms.get()]
        while(new_idx != []):
            atom_idx = new_idx.pop(0)
            if atom_idx not in atom_list_idx:
                atom_list_idx.append(atom_idx)
            connected_atoms = list(np.where(graph[atom_idx])[0])
            for connected_atom_idx in connected_atoms:
                if connected_atom_idx not in atom_list_idx:
                    if atom_idx in cen_atoms_idx and connected_atom_idx in mea_atoms_idx:
                        new_idx.append(connected_atom_idx)
                    elif atom_idx in mea_atoms_idx and connected_atom_idx in cen_atoms_idx:
                        new_idx.append(connected_atom_idx)
        return atom_list_idx
