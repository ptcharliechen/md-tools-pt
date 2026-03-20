from kit.args import args
arg = args({"output": "split", "incar": [None, str, "INCAR path"], "kpoints": [None, str, "KPOINTS path"], "potcar": [None, str, "POTCAR path"]}, \
        "POSCAR", steps=False, output_isdir=True)
from os import path, mkdir
from shutil import copy
from copy import deepcopy
from kit.fundamental import Atom
from kit.vasp import POSCAR
from kit.interface import line

if __name__ == "__main__":
    ads_poscar = POSCAR(arg)
    print("\nAdsorbate:")
    ads_atoms = line(ads_poscar)
    ads_poscar.read_all(ads_atoms)
    while(1):
        inp = input("\nDo you want to use the remaining as bulk atoms (y/n) [y]: ")
        if inp == "":
            inp = 'y'
        if inp.lower() == 'n':
            print("Bulk:")
            bulk_poscar = POSCAR(arg)
            bulk_atoms = line(bulk_poscar)
            bulk_poscar.read_all(bulk_atoms)
            break
        elif inp.lower() == 'y':
            bulk_atoms = deepcopy(Atom(ads_poscar))
            ads_atom_list = ads_atoms.get()
            for atom in range(len(ads_poscar.elements)):
                if atom not in ads_atom_list:
                    bulk_atoms.put(atom)
            bulk_poscar = POSCAR(arg)
            bulk_poscar.read_all(bulk_atoms)
            break
    total_poscar = POSCAR(arg)
    total_atoms = deepcopy(Atom(arg))
    total_atoms.put(ads_atoms.get())
    total_atoms.put(bulk_atoms.get())
    total_atoms.sort()
    total_poscar.read_all(total_atoms)
    for kind in ["incar", "kpoints", "potcar"]:
        if getattr(arg.args, kind) is None:
            continue
        if path.isfile(path.join(getattr(arg.args, kind), kind.upper())):
            setattr(arg.args, kind, path.join(getattr(arg.args, kind), kind.upper()))
        elif not path.isfile(getattr(arg.args, kind)):
            print("Warning: {} is invalid.".format(kind))
            setattr(arg.args, kind, None)
    if (arg.args.incar is None) or (arg.args.kpoints is None) or (arg.args.potcar is None):
        while(1):
            yn = input("\nDo you want to take INCAR, KPOINTS & POTCAR in directories? (y/n) [y]: ")
            if yn == "":
                yn = 'y'
            if yn.lower() == 'y':
                for kind in ["incar", "kpoints", "potcar"]:
                    if getattr(arg.args, kind) is None:
                        pass
                    elif path.isfile(getattr(arg.args, kind)):
                        continue
                    while(1):
                        tmp = input("{} path: ".format(kind.upper()))
                        if path.isfile(tmp):
                            setattr(arg.args, kind, path.abspath(tmp))
                            break
                        elif path.isfile(path.join(getattr(arg.args, kind), kind.upper())):
                            setattr(arg.args, kind, path.join(path.abspath(tmp), kind.upper()))
                            break
                        else:
                            print("This file doesn't exist.\n")
                break
            elif yn.lower() == 'n':
                break
            else:
                print("Warning: Input 'y' or 'n'.")
    else:
        yn = 'y'

    mkdir(arg.args.output)
    elements = []
    for poscar, atoms, kind in zip([bulk_poscar, ads_poscar, total_poscar], [bulk_atoms, ads_atoms, total_atoms], ["bulk", "adsorbate", "total"]):
        mkdir(path.join(arg.args.output, kind))
        write_poscar = POSCAR(arg)
        write_poscar.args.output = path.join(arg.args.output, kind, "POSCAR")
        write_poscar.lattice = poscar.lattice
        element = []
        tmp = atoms.get()
        for idx, _ in enumerate(poscar.elements):
            if idx in tmp:
                element.append(_)
        elements.append(element)
        write_poscar.elements = element
        write_poscar.atoms = atoms
        write_poscar.fractional_position = poscar.fractional_position
        write_poscar.write_all()

        if yn.lower() == 'y':
            copy(arg.args.incar, path.join(arg.args.output, kind))
            copy(arg.args.kpoints, path.join(arg.args.output, kind))
    
    if yn.lower() == 'y':
        with open(path.join(arg.args.output, "total", "POTCAR"), "w") as store_file:
            flag = False
            with open(arg.args.potcar, "r") as open_file:
                for line in open_file:
                    if not flag and "PAW" in line:
                        if line.split()[1].split("_")[0] in elements[2]:
                            flag = True
                    if "End of Dataset" in line and flag:
                        store_file.write(line)
                        flag = False
                        continue
                    elif flag:
                        store_file.write(line)

        with open(path.join(arg.args.output, "bulk", "POTCAR"), "w") as store_file:
            flag = False
            with open(arg.args.potcar, "r") as open_file:
                for line in open_file:
                    if not flag and "PAW" in line:
                        if line.split()[1].split("_")[0] in elements[0]:
                            flag = True
                    if "End of Dataset" in line and flag:
                        store_file.write(line)
                        flag = False
                    elif flag:
                        store_file.write(line)

        with open(path.join(arg.args.output, "adsorbate", "POTCAR"), "w") as store_file:
            flag = False
            with open(arg.args.potcar, "r") as open_file:
                for line in open_file:
                    if not flag and "PAW" in line:
                        if line.split()[1].split("_")[0] in elements[1]:
                            flag = True
                    if "End of Dataset" in line and flag:
                        store_file.write(line)
                        flag = False
                    elif flag:
                        store_file.write(line)
    print("Done!")
