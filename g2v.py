from kit.args import args
arg = args({"output": "POSCAR"}, "POSCAR", steps=False, atoms=False)
import re, warnings
import numpy as np
from kit.software import Gaussian
from kit.vasp import POSCAR
from kit.args import args
warnings.filterwarnings("ignore", category=DeprecationWarning)

if __name__ == "__main__":
    vector = []
    gaussian = Gaussian(arg)
    gaussian.read_all()
    if len(gaussian.lattice) == 0:
        flag = False
        x, y, z = gaussian.cartesian_position[:, 0], gaussian.cartesian_position[:, 1], gaussian.cartesian_position[:, 2]
        while(1):
            size = input("Lattice size (s/m/l/u): ")
            if size.lower() == "s":
                factor = 1
                break
            elif size.lower() == "m":
                factor = 2
                break
            elif size.lower() == "l":
                factor = 6
                break
            elif size.lower() == "u":
                from re import compile
                tmp = input("Lattice shape: ")
                if compile(r"^(\d+(\.\d+)?)(,(\d+(\.\d+)?)){2}(\s(\d+(\.\d+)?)(,(\d+(\.\d+)?)){2}){2}$").match(tmp) is None or len(tmp.split()) != 3:
                    print("Warning: Three numbers separated by cammas as a set. Three sets separated by spaces, like '3,0,0 0,3,0 0,0,3'.")
                    continue
                else:
                    tmp = tmp.split()
                    for idx, val in enumerate(tmp):
                        vector.append([float(num) for num in val.split(",")])
                    flag = True
                    break
            else:
                print("Warning: Input error.")
        max_dist = 0
        for i in range(len(x)):
            for j in range(len(y)):
                dist = np.linalg.norm(np.array([x[i], y[i], z[i]]) - np.array([x[j], y[j], z[j]]))
                if dist > max_dist:
                    max_dist = dist
        x += abs(min(x))
        y += abs(min(y))
        z += abs(min(z))
        if not flag:
            x += max_dist*factor
            y += max_dist*factor
            z += max_dist*factor
            if np.mean(np.array(z))*2 < max([np.mean(np.array(x))*2, np.mean(np.array(y))*2])/2.5:
                z = (x+y)/2
            x_mean = np.mean(np.array(x))*2
            y_mean = np.mean(np.array(y))*2
            z_mean = np.mean(np.array(z))*2
            max_value = max([x_mean, y_mean, z_mean])
            vector = np.array([[max_value, 0, 0], [0, max_value, 0], [0, 0, max_value]])
        else:
            x += vector[0][0]/2 + vector[1][0]/2 + vector[2][0]/2
            y += vector[0][1]/2 + vector[1][1]/2 + vector[2][1]/2
            z += vector[0][2]/2 + vector[1][2]/2 + vector[2][2]/2
        gaussian.lattice = vector
        pos = np.zeros(gaussian.cartesian_position.shape)
        pos[:, 0], pos[:, 1], pos[:, 2] = x, y, z
        gaussian.cartesian_position = pos
    poscar = POSCAR(arg)
    poscar.bridge = gaussian
    while(1):
        tmp = input("Fixed atom border [-1]: ")
        if tmp == "":
            tmp = "-1"
        if re.compile(r"^-?\d+\.?\d*$").match(tmp) is None:
            print("Warning: Input a decimal.")
        elif float(tmp) > 1 or float(tmp) < -1:
            print("Warning: Input a number between -1 and 1.")
        else:
            poscar.relax = True
            poscar.relax_border = float(tmp)
            break

    poscar.write_all()
    print("Done!")
