from os import path
from sys import argv
from re import compile
from numpy import zeros
from kit.fundamental import Lattice
from kit.vasp import POSCAR
from kit.interface import args

if __name__ == "__main__":
    arg = args({"output": "POSCAR"}, "data.lmp")
    if len(argv) < 3 and path.isfile("elements"):
        argv.append("elements")
        with open("elements") as read_file:
            tmp = read_file.read().replace('\n', '')
            if compile(r"^[a-zA-Z]+( [a-zA-Z]+)*$").match(tmp) is not None:
                argv.append(tmp)
                elementsData = tmp.split()
    if len(argv) < 3:
        while(len(argv) < 3):
            filename = input("Input elements file: ")
            if path.isfile(filename):
                with open(filename) as read_file:
                    tmp = read_file.read().replace('\n', '')
                    if compile(r"^[a-zA-Z]+( [a-zA-Z]+)*$").match(tmp) is not None:
                        argv.append(filename)
                        elementsData = tmp.split()
            else:
                print("Warning: File not found.")
    
    AtomInfo = {}
    with open(argv[1]) as read_file:
        for line in read_file:
            if "atoms" in line:
                atom_sum = int(line.split()[0])
            elif "xlo" in line:
                lattice = Lattice()
                i = 0
                tmp = zeros((3, 3))
                sp = line.split()
                CellLowerX, CellUpperX = float(sp[0]), float(sp[1])
                tmp[i, i] = CellUpperX - CellLowerX
                sp = next(read_file).split()
                i += 1
                CellLowerY, CellUpperY = float(sp[0]), float(sp[1])
                tmp[i, i] = CellUpperY - CellLowerY
                sp = next(read_file).split()
                i += 1
                CellLowerZ, CellUpperZ = float(sp[0]), float(sp[1])
                tmp[i, i] = CellUpperZ - CellLowerZ
                lattice.lattice = tmp
            elif "Atoms" in line:
                line = next(read_file)
                for i in range(atom_sum):
                    sp = next(read_file).split()
                    AtomInfo[int(sp[0])-1] = [elementsData[int(sp[2])-1], float(sp[4])+CellLowerX, float(sp[5])+CellLowerY, float(sp[6])+CellLowerZ]
    positions, elements, num = [], [], []
    for i in sorted(AtomInfo.keys()):
        num.append(i)
        elements.append(AtomInfo[i][0])
        positions.append(AtomInfo[i][1:])
    
    poscar = POSCAR(arg)
    poscar.lattice = lattice.lattice
    poscar.rearrange(elements, positions, num, typ='c')
    poscar.write_all(typ='f')
    
    print("Done!")
