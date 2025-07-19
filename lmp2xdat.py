import sys
from numpy import zeros, array
from kit.fundamental import Lattice
from kit.accelerate import c2f

def rearange(AtomInfo):
    IndexDict, EleNum = {}, {}
    elements = []
    
    for val in AtomInfo.values():
        if val[0] not in elements:
            elements.append(val[0])
    idx = 0
    for element in elements:
        for key, val in AtomInfo.items():
            if element == val[0]:
                IndexDict[idx] = key
                EleNum[element] = 1 if element not in EleNum.keys() else EleNum[element] + 1
                idx += 1
    return EleNum, IndexDict

if __name__ == '__main__':
    lattice = Lattice()
    with open(sys.argv[1]) as read_file:
        FirstFlag = True
        AtomInfo = {}
        step = 1
        with open("XDATCAR", "w") as write_file:
            for line in read_file:
                if "ITEM: TIMESTEP" in line:
                    line = next(read_file)
                    line = next(read_file)
                    line = next(read_file)
                    AtomSum = int(line.split()[0])
                    line = next(read_file)
                    sp = next(read_file).split()
                    tmp = zeros((3, 3))
                    i = 0
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
                    line = next(read_file)
                    for i in range(AtomSum):
                        sp = next(read_file).split()
                        AtomInfo[int(sp[0])] = [sp[3], float(sp[5])+CellLowerX, float(sp[6])+CellLowerY, float(sp[7])+CellLowerZ]
                    if FirstFlag or lattice.NpT_flag:
                        write_file.write("trajectory\n  1.0\n")
                        lat = lattice.lattice if lattice.lattice.ndim == 2 else lattice.lattice[-1]
                        for l in lat:
                            write_file.write(f"   {l[0]:.10f}   {l[1]:.10f}   {l[2]:.10f}\n")
                        if FirstFlag:
                            FirstFlag = False
                            EleNum, IndexDict = rearange(AtomInfo)
                        write_file.write("".join(f"   {idx}" for idx in EleNum.keys()))
                        write_file.write('\n')
                        write_file.write("".join(f"   {idx}" for idx in EleNum.values()))
                        write_file.write('\n')
                    write_file.write(f"Direct configuration=    {step}\n")
                    positions = []
                    for i in range(AtomSum):
                        positions.append(array(AtomInfo[IndexDict[i]][1:]))
                    positions = c2f(array(positions), lat)
                    for i in range(AtomSum):
                        write_file.write(f"   {positions[i][0]:.8f}    {positions[i][1]:.8f}    {positions[i][2]:.8f}\n")
                    step += 1
    print("Done!")
