from os import path, getcwd, system, listdir
from sys import argv
from numpy import array_split
from kit.machine import Slurm_VASP

python_script = """
from getpass import getuser
from os import remove
import sys
sys.path.append(f"/home/{getuser()}/script/kit")
from machine import VASP

multiprocess = VASP(sub)
multiprocess.processes = process_num
multiprocess.build_panel(dir)
multiprocess.parallel(processes=num)
"""

if __name__ == "__main__":
    while(1):
        inp = input("Do you want to submit all jobs? (y/n) [y]: ")
        if inp == '':
            inp = 'y'
        if inp.lower() == 'y':
            dirs = []
            for dir in listdir(getcwd()):
                if path.isdir(path.join(getcwd(), dir)):
                    for file in ["INCAR", "POSCAR", "POTCAR", "KPOINTS"]:
                        if not path.isfile(path.join(getcwd(), dir, file)):
                            print(f"Warning: '{file}' can not be found in {dir}.")
                            break
                    else:
                        dirs.append(dir)
            if dirs == []:
                print("Error: No jobs to submit.")
                exit(1)
            else:
                dirs = sorted(dirs)
            break
        elif inp.lower() == 'n':
            dirs = []
            idx = 1
            while(1):
                print("Input directory name ('end' to finish):")
                dir = input(f"{idx}: ")
                if dir.lower() == "end":
                    break
                elif path.isdir(path.join(getcwd(), dir)):
                    for file in ["INCAR", "POSCAR", "POTCAR", "KPOINTS"]:
                        if not path.isfile(path.join(getcwd(), dir, file)):
                            print(f"Warning: '{file}' can not be found in {dir}.")
                            break
                    else:
                        dirs.append(dir)
                        idx += 1
                else:
                    print("Warning: Directory not found.")
            break
        else:
            print("Warning: Input error.")

    slurm = Slurm_VASP()
    slurm.name = name = argv[1] if len(argv) > 1 else "submit"
    slurm.input_cores()
    slurm.input_process(default=1)
    slurm.input_nodes()
    slurm.input_partition()
    
    dirs = array_split(dirs, slurm.nodes)
    for idx, dir_split in enumerate(dirs):
        with open(path.join(getcwd(), f"submit_{idx}.py"), "w") as write_file:
            write_file.write(python_script.replace("process_num", str(int(slurm.cores/slurm.processes))).replace("sub", str(list(dir_split))).replace("dir", f"\".panel_{idx}.csv\"").replace("num", str(slurm.processes))+f"remove(\"submit_{idx}.py\")\n")
        system(f"chmod 700 submit_{idx}.py")
    
    slurm_script = []
    for idx, line in enumerate(slurm.slurm.split("\n")):
        if ".log" in line:
            command = [f"srun -N1 -n1 -w ${{nodes[{node}]}} {getcwd()}/submit_{node}.py &" for node in range(slurm.nodes)]
            slurm_script.append("\n".join(command) + '\n' + f"wait >& {name}.log")
        elif "$SLURM_JOB_NODELIST" in line:
            slurm_script.append(line)
            slurm_script.append("nodes=($(scontrol show hostnames $SLURM_NODELIST))")
        elif "vdw_kernel.bindat" in line:
            continue
        else:
            slurm_script.append(line)
    
    for node in range(slurm.nodes):
        slurm_script.append(f"rm .panel_{node}.csv")
    slurm.slurm = slurm_script
    slurm.write_slurm()

    system(f"chmod 755 {name}.slurm")
    system(f"sbatch {name}.slurm")
