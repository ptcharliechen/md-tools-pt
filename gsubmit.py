from os import path, getcwd, system, listdir
from numpy import array_split
from kit.machine import Partition

if __name__ == "__main__":
    while(1):
        inp = input("Do you want to submit all jobs? (y/n) [y]: ")
        if inp == '':
            inp = 'y'
        if inp.lower() == 'y':
            files = []
            for file in listdir(getcwd()):
                if file.endswith(".gjf"):
                    files.append(file)
            if files == []:
                print("Error: No jobs to submit.")
                exit(1)
            break
        elif inp.lower() == 'n':
            files = []
            idx = 1
            print("Input file name ('end' to finish):")
            while(1):
                file = input(f"{idx}: ")
                if file.lower() == "end":
                    break
                elif not path.isfile(path.join(getcwd(), file)):
                    print("Warning: File not found.")
                elif file in files:
                    print("Warning: File in list.")
                elif file.endswith(".gjf"):
                    files.append(file)
                    idx += 1
                else:
                    print("Warning: Import gjf format.")
            break
        else:
            print("Warning: Input error.")

    partition = Partition()
    partition.input_resource(process=True)
    
    for file in files:
        system(f"/home/j14jjc00/bin/gpreprocess.py {path.join(getcwd(), file)} {int(partition.cores/partition.processes)}")

    python_script = f"""from os import remove
from machine import Gaussian
multiprocess = Gaussian(sub)
multiprocess.build_panel(file)
multiprocess.parallel(processes={partition.processes})
"""
    
    files = array_split(sorted(files), partition.nodes)
    for idx, file_split in enumerate(files):
        with open(path.join(getcwd(), f"submit_{idx}.py"), "w") as write_file:
            write_file.write(python_script.replace("sub", str(list(file_split))).replace("file", f"\".panel_{idx}.csv\"")+f"remove(\"submit_{idx}.py\")\n")
        system(f"chmod 700 submit_{idx}.py")
    
    slurm_script = f"""#!/bin/bash
#SBATCH -J submit
#SBATCH -o submit.out
#SBATCH -e submit.err
#SBATCH -p {partition.partition}
#SBATCH --nodes {partition.nodes}
#SBATCH --ntasks-per-node={partition.cores}
#SBATCH --export=all
#SBATCH --no-requeue

cd $SLURM_SUBMIT_DIR
echo $SLURM_JOB_NODELIST > NODE.$SLURM_JOBID
nodes=($(scontrol show hostnames $SLURM_NODELIST))

module purge
# Gaussian path
echo -n \"g16  : \" && which g16

ulimit -s unlimited
sleep 0.5
echo \"Begin time: `date`\"
script >& submit.log
echo \"Finish time: `date`\"
"""
    
    command, remove_command = "", ""
    for node in range(partition.nodes):
        command += f"srun -N1 -n1 -w ${{nodes[{node}]}} {getcwd()}/submit_{node}.py &\n"
        remove_command += f"rm .panel_{node}.csv\n"
    command += "wait"
    slurm_script = slurm_script.replace("script", command)
    with open(path.join(getcwd(), "submit.slurm"), "w") as slurm_file:
        slurm_file.write(slurm_script+remove_command)
    system("chmod 755 submit.slurm")
    system("sbatch submit.slurm")
