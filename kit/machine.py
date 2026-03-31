from os import path, system, cpu_count, popen
from abc import ABC, abstractmethod
from glob import glob

class Partition:
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._partitions = ["developement", "ct112", "ct448", "ct1k", "ct2k", "ct4k", "ct8k"]
        self._cores, self._processes, self._nodes = cpu_count(), 1, 1
        self._partition = "developement"
    @property
    def cores(self):
        return self._cores
    @property
    def processes(self):
        return self._processes
    @property
    def nodes(self):
        return self._nodes
    @property
    def partition(self):
        return self._partition
    @partition.setter
    def partition(self, partition):
       if partition not in self._partitions:
           raise Exception(f"The partition '{partition}' is not available.")
       self._partition = partition
    @cores.setter
    def cores(self, cores):
        if not str(cores).isdigit():
            raise Exception("The format should be a positive integer.")
        elif not isinstance(cores, int):
            self._cores = int(cores)
        else:
            self._cores = cores
    @processes.setter
    def processes(self, processes):
        if not str(processes).isdigit():
            raise Exception("The format should be a positive integer.")
        elif not isinstance(processes, int) and int(processes) <= self._cores:
            self._processes = int(processes)
        elif isinstance(processes, int) and processes <= self._cores:
            self._processes = processes
        else:
            raise ValueError(f"The number of processes should be smaller than {self._cores}.")
    @nodes.setter
    def nodes(self, nodes):
        if not str(nodes).isdigit():
            raise Exception("The format should be a positive integer.")
        elif not isinstance(nodes, int):
            self._nodes = int(nodes)
        else:
            self._nodes = nodes
    def input_cores(self, default=cpu_count()):
        while(1):
            tmp = input(f"Number of cores per node [{default}]: ")
            if tmp == "":
                tmp = default
            if str(tmp).isdigit() and int(tmp) <= cpu_count():
                self._cores = int(tmp)
                break
            elif str(tmp).isdigit() and int(tmp) > cpu_count():
                print("Warning: Input a positive integer below {cpu_cpunt()")
            else:
                print("Warning: The format should be a positive integer.")
    def input_nodes(self, default=1):
        while(1):
            tmp = input(f"Number of nodes [{default}]: ")
            if tmp == "":
                tmp = default
            if str(tmp).isdigit() and self._cores * int(tmp) > 8400:
                print(f"Warning: Input a positive integer below {tmp}.")
            elif str(tmp).isdigit():
                self._nodes = int(tmp)
                break
            else:
                print("Warning: The format should be a positive integer.")
    def input_process(self, default=-1):
        if default < 1:
            default = int(self._cores / 8)
        while(1):
            tmp = input(f"Number of processes per node [{default}]: ")
            if tmp == "":
                tmp = default
            if str(tmp).isdigit() and int(tmp) <= self._cores:
                self._processes = int(tmp)
                break
            elif str(tmp).isdigit() and int(tmp) > self._cores:
                print(f"Warning: Input a positive integer below {self._cores}")
            else:
                print("Warning: The format should be a positive integer.")
    def input_partition(self):
        while(1):
            tmp = input("'development' partition or not (y/n) [n]: ")
            if tmp == "":
                tmp = "n"
            if tmp.lower() == "y":
                self._partition = "development"
                break
            elif tmp.lower() == "n":
                if self._cores * self._nodes < 113:
                    self._partition = "ct112"
                elif 112 < self._cores * self._nodes < 449:
                    self._partition = "ct224"
                elif 448 < self._cores * self._nodes < 1121:
                    self._partition = "ct1k"
                elif 1120 < self._cores * self._nodes < 2241:
                    self._partition = "ct2k"
                elif 2240 < self._cores * self._nodes < 4481:
                    self._partition = "ct4k"
                elif 4480 < self._cores * self._nodes < 8961:
                    self._partition = "ct8k"
                break
            else:
                print("Warning: Input 'y' or 'n'.")
    def input_resource(self, process=False):
        self.input_cores()
        if process:
            self.input_process()
        self.input_nodes()
        self.input_partition()

class Job:
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        from os import getcwd
        self._format = "vasp641"
        self._name = "default"
        self._dirpath = getcwd()
    @property
    def format(self):
        return self._format
    @property
    def dir_path(self):
        return self._dirpath
    @property
    def name(self):
        return self._name
    @format.setter
    def format(self, format):
        self._format = format
    @dir_path.setter
    def dir_path(self, jobpath):
        if path.isfile(jobpath):
            self._dirpath = path.abspath(path.dirname(jobpath))
            self._name = path.basename(jobpath).split(".")[0]
        elif path.isdir(jobpath):
            self._dirpath = path.abspath(jobpath)
        else:
            raise FileNotFoundError(f"{jobpath} does not exist.")
    @name.setter
    def name(self, name):
        self._name = name.split(".")[0]
    def input_name(self, same_name=True):
        tmp = input("\tJob name: ")
        if tmp == "":
            tmp = "default"
        i = 1
        while(same_name and path.isfile(path.join(self._dirpath, tmp))):
            tmp = tmp.split("_")[0] + "_" + str(i)
            i += 1
        self._name = tmp.split(".")[0]

class Slurm(Job):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._slurm = ""
    @property
    def slurm(self):
        return self._slurm
    @slurm.setter
    def slurm(self, slurm):
        if isinstance(slurm, str):
            self._slurm = slurm
        elif isinstance(slurm, list) or isinstance(slurm, tuple):
            self._slurm = "".join(line if line.endswith('\n') else line + '\n' for line in slurm)
        else:
            raise TypeError(f"{type(slurm)} is not supported.")
    def write_slurm(self):
        with open(path.join(self._dirpath, self._name+".slurm"), "w") as write_file:
            write_file.write(self.slurm)

class Slurm_VASP(Slurm, Partition):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
    @Slurm.slurm.getter
    def slurm(self):
        if self._slurm != "":
            return super().slurm
        elif self._format == "vasp641":
            self.vasp641()
        elif self._format == "vasp611":
            self.vasp611()
        elif self._format == "vasp544":
            self.vasp544()
        else:
            raise NotImplementedError("'vasp' format is applied only.")
        return self._slurm
    def vasp641(self):
        self._slurm = f"""#!/bin/bash

### Executed by VASP.6.4.1

#SBATCH -J {self._name}
#SBATCH -o {self._name}.out
#SBATCH -e {self._name}.err
#SBATCH -p {self._partition}
#SBATCH --nodes={self._nodes}
#SBATCH --ntasks-per-node={self._cores}
#SBATCH --export=all
#SBATCH --no-requeue

cd $SLURM_SUBMIT_DIR
echo $SLURM_JOB_NODELIST > NODE.$SLURM_JOBID

# source intel compiler
module purge
# intel path

# source vasp
# VASP path

# var check
echo -n "Intel mpiifort : " && which mpiifort
echo -n "Intel MPI      : " && which mpiexec.hydra
echo -n "vasp_std       : " && which vasp_std

ulimit -s unlimited
unset I_MPI_PMI_LIBRARY 
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0
sleep 2

echo "Begin time  : `date`"
mpiexec.hydra -np $SLURM_NTASKS -ppn $SLURM_NTASKS_PER_NODE vasp_std >& {self._name}.log
echo "Finish time : `date`"
"""
    def vasp611(self):
        self._slurm = f"""#!/bin/bash

### submit vasp.6.1.1

#SBATCH -J {self._name}
#SBATCH -o {self._name}.out
#SBATCH -e {self._name}.err
#SBATCH -p {self._partition}
#SBATCH --nodes={self._nodes}
#SBATCH --ntasks-per-node={self._cores}
#SBATCH --export=all
#SBATCH --no-requeue

cd $SLURM_SUBMIT_DIR
echo $SLURM_JOB_NODELIST > NODE.$SLURM_JOBID

# source intel compiler
module purge
# intel path

# source vasp
# VASP path

# var check
echo -n "Intel mpiifort : " && which mpiifort
echo -n "Intel MPI      : " && which mpiexec.hydra
echo -n "vasp_std       : " && which vasp_std

ulimit -s unlimited
unset I_MPI_PMI_LIBRARY 
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0
sleep 2

echo "Begin time  : `date`"
mpiexec.hydra -np $SLURM_NTASKS -ppn $SLURM_NTASKS_PER_NODE vasp_std >& {self._name}.log
echo "Finish time : `date`"
"""
    def vasp544(self):
        self._slurm = f"""#!/bin/bash

### submit vasp.5.4.4

#SBATCH -J {self._name}
#SBATCH -o {self._name}.out
#SBATCH -e {self._name}.err
#SBATCH -p {self._partition}
#SBATCH --nodes={self._nodes}
#SBATCH --ntasks-per-node={self._cores}
#SBATCH --export=all
#SBATCH --no-requeue

cd $SLURM_SUBMIT_DIR
echo $SLURM_JOB_NODELIST > NODE.$SLURM_JOBID

# source intel compiler
module purge
# intel path

# source vasp
# VASP path

# var check
echo -n "Intel mpiifort : " && which mpiifort
echo -n "Intel MPI      : " && which mpiexec.hydra
echo -n "vasp_std       : " && which vasp_std

ulimit -s unlimited
unset I_MPI_PMI_LIBRARY 
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0
sleep 2

echo "Begin time  : `date`"
mpiexec.hydra -np $SLURM_NTASKS -ppn $SLURM_NTASKS_PER_NODE vasp_std >& {self._name}.log
echo "Finish time : `date`"
"""

class Slurm_QE(Slurm, Partition):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
    def qe6(self):
        if path.isfile(path.join(self._dirpath, "QE.in")):
            FILE = "QE.in"
            SCRIPT = "pw.x"
        elif path.isfile(path.join(self._dirpath, "RISM.in")):
            FILE = "RISM.in"
            SCRIPT = "pprism.x"
        self._slurm = f"""#!/bin/bash

### Executed by QE v6.7

#SBATCH -J {self._name}
#SBATCH -o {self._name}.out
#SBATCH -e {self._name}.err
#SBATCH -p {self._partition}
#SBATCH --nodes={self._nodes}
#SBATCH --ntasks-per-node={self._cores}
#SBATCH --export=all
#SBATCH --no-requeue

cd $SLURM_SUBMIT_DIR
echo $SLURM_JOB_NODELIST > NODE.$SLURM_JOBID

# source intel compiler
module purge
# intel path

# source QE v.6-7
# QE path

# var check
echo -n "Intel mpiifort : " && which mpiifort
echo -n "Intel MPI      : " && which mpiexec.hydra
echo -n "QE             : " && which {SCRIPT}

ulimit -s unlimited
unset I_MPI_PMI_LIBRARY 
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0
sleep 2

echo "Begin time  : `date`"
mpiexec.hydra -np $SLURM_NTASKS -ppn $SLURM_NTASKS_PER_NODE {SCRIPT} < {FILE} >& {self._name}.log
echo "Finish time : `date`"
"""
    def qe7(self):
        if path.isfile(path.join(self._dirpath, "QE.in")):
            FILE = "QE.in"
            SCRIPT = "pw.x"
        elif path.isfile(path.join(self._dirpath, "RISM.in")):
            FILE = "RISM.in"
            SCRIPT = "pprism.x"
        self._slurm = f"""#!/bin/bash

### Executed by QE v7.1

#SBATCH -J {self._name}
#SBATCH -o {self._name}.out
#SBATCH -e {self._name}.err
#SBATCH -p {self._partition}
#SBATCH --nodes={self._nodes}
#SBATCH --ntasks-per-node={self._cores}
#SBATCH --export=all
#SBATCH --no-requeue

cd $SLURM_SUBMIT_DIR
echo $SLURM_JOB_NODELIST > NODE.$SLURM_JOBID

# source intel compiler
module purge
# intel path

# source QE v.7-1
# QE path

# var check
echo -n "Intel mpiifort : " && which mpiifort
echo -n "Intel MPI      : " && which mpiexec.hydra
echo -n "QE             : " && which {SCRIPT}

ulimit -s unlimited
unset I_MPI_PMI_LIBRARY 
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0
sleep 2

echo "Begin time  : `date`"
mpiexec.hydra -np $SLURM_NTASKS -ppn $SLURM_NTASKS_PER_NODE {SCRIPT} < {FILE} >& {self._name}.log
echo "Finish time : `date`"
"""
    @Slurm.slurm.getter
    def slurm(self):
        if self._slurm != "":
            return super().slurm
        elif self._format == "qe6":
            self.qe6()
        elif self._format == "qe7":
            self.qe7()
        else:
            raise NotImplementedError("'qe' format is applied only.")
        return self._slurm

class Slurm_Conquest(Slurm, Partition):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._format = "conquest"
    def conquest(self):
        total_processes = self._processes * self._nodes
        if path.isfile(path.join(self._dirpath, "ionpos.dat")):
            line = len(open(path.join(self._dirpath, "ionpos.dat"),'r').readlines()) - 4
        else:
            raise FileNotFoundError(f"Error: The 'ionpos.dat' file is not found in {self._dirpath}.")
        if total_processes > line:
            raise Exception(f"Error: The number of used processes should be smaller than the number of atoms ({total_processes} processes > {line} atoms). Decrease the process number.")
        self._slurm = f"""#!/bin/bash

### Executed by Conquest

#SBATCH -J {self._name}
#SBATCH -o {self._name}.out
#SBATCH -e {self._name}.err
#SBATCH -p {self._partition}
#SBATCH --nodes={self._nodes}
#SBATCH --ntasks-per-node={self._cores}
#SBATCH --export=all
#SBATCH --no-requeue

cd $SLURM_SUBMIT_DIR
echo $SLURM_JOB_NODELIST > NODE.$SLURM_JOBID

# source intel compiler
module purge
# intel path

# source Conquest
# Conquest path

# var check
echo -n "intel mpiifort : " && which mpiifort
echo -n "intel mpi      : " && which mpiexec.hydra
echo -n "Conquest       : " && which Conquest

ulimit -s unlimited
unset I_MPI_PMI_LIBRARY 
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0
sleep 2

echo "Begin time  : `date`"
mpiexec.hydra -np {total_processes} -ppn {self._processes} Conquest < Conquest_input >& {self._name}.log
echo "Finish time : `date`"
"""
    @Slurm.slurm.getter
    def slurm(self):
        if self._slurm != "":
            return super().slurm
        else:
            self.conquest()
            return self._slurm

def input_format():
    while(1):
        tmp = input("\tType (v64: VASP641, v61: VASP611, v5: VASP544, c: Conquest, q6: QE6, q7: QE7) [v64]: ")
        if tmp == "" or tmp.lower() == "v64":
            return "vasp641"
        elif tmp.lower() == "v61":
            return "vasp611"
        elif tmp.lower() == "v5":
            return "vasp544"
        elif tmp.lower() == "c":
            return "conquest"
        elif tmp.lower() == "q6":
            return "qe6"
        elif tmp.lower() == "q7":
            return "qe7"
        else:
            print("Warning: Input 'v64' or 'v61' or 'v5' or 'c' or 'q6' or 'q7'.")

def submit_slurm(SlurmPath, SlurmFile):
    from os import popen
    with popen(f"cd {SlurmPath} && sbatch {SlurmFile}") as read_command:
        content = read_command.read()
        if "Submitted batch job" in content:
            return content.split()[3], content.replace('\n', '')
        else:
            return "-1", content

def parallel(scripts, *args, cpu=int(cpu_count()/2.5), compiler="python3") -> None:
    """
    Refactored controller function to use ThreadPoolExecutor for parallel execution.
    
    Args:
        self: the instance of the class
        max_jobs (int): maximum number of jobs for ThreadPoolExecutor
    
    Returns:
        None
    """
    import subprocess
    from concurrent.futures import ThreadPoolExecutor, as_completed
    from os import getcwd, remove
    cmds = []
    for i, script, arg in zip(tuple(range(len(scripts))), scripts, args[0]):
        with open(f"{i}.py", "w") as write_file:
            write_file.write(script)
        cmd = f"{compiler} {path.join(getcwd(), f'{i}.py')}"
        system(f"chmod 700 {path.join(getcwd(), f'{i}.py')}")
        for parameter in arg:
            cmd += f" {parameter}"
        cmds.append(cmd)
    
    # Create ThreadPoolExecutor with specified number of workers
    with ThreadPoolExecutor(cpu) as executor:
        # Dictionary to store futures and corresponding job IDs
        futures = {}
        
        # Construct command for execution
        for cmd in cmds:
            # Submit command to ThreadPoolExecutor and store future with job ID
            futures = [executor.submit(lambda cmd: subprocess.run(cmd, shell=True), cmd)]

        # Iterate over completed futures
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as err:
                print(f"System error: {err}")
    
    for i in range(len(scripts)):
        remove(path.join(getcwd(), f"{i}.py"))

class Multiprocess(ABC):
    def __init__(self, file_list: list, **kwargs) -> None:
        """
        Initialize Multiprocess with a file list.
        
        Args:
            file_list (list): List of files to be processed.
        """
        super().__init__(**kwargs)
        from os import getcwd
        self._format = None
        self._dir_path = getcwd()
        # Initialize job dictionary
        self._job_dict = {}
        
        # Check and create directories
        self._dir_check()

        # Store the file list
        self._file_list = file_list
        self._cmds = None
    @abstractmethod
    def _move_to_queue(self, file_list):
        """
        Move the jobs to 'queue' directory.
        """
        for file in file_list:
            system(f"mv {path.join(self._dir_path, file)} {path.join(self._dir_path, 'queue')}")
    @abstractmethod
    def _default_cmd(self):
        """
        If self._cmds is None, this function will be called as deafult commands.
        """
        pass
    @abstractmethod
    def _exception(self, job_id, e):
        """
        Handle the exception if the process has crashed.
        Stamp status, write the error in 'error message' to the panel, and move the job to 'error'.
        """
        self._file_stat[job_id]["status"] = "Error"
        base_file = path.join(self._dir_path, "queue", self._file_stat[job_id]["file"])
        system(f"mv {base_file}.gjf {path.join(self._dir_path, 'error')}/")
        if path.isfile(base_file + ".log"):
            system(f"mv {base_file}.log {path.join(self._dir_path, 'error')}/")
        self._file_stat[job_id]["error message"] = str(exception)
    @abstractmethod
    def _end(self, job_id):
        """
        After the process is done, operate this function.
        Even if the process is done without any errors, the job could give some errors. Detect these errors in log files, write them in 'wrong message' to the panel, and send the job to 'wrong'.
        If the job completes normally, send the job to 'done'.
        """
        pass
    @property
    def commands(self):
        return self._cmds
    @commands.setter
    def commands(self, commands):
        assert len(commands) == len(self._file_list), "The number of commands does not match the number of files."
        self._cmds = commands
    def stamp_and_run(self, cmd: str, job_id: int) -> None:
        """
        Executes a command and measures its execution time.

        Args:
            cmd (str): The command to be executed.

        Returns:
            None
        """
        from time import strftime, localtime, time

        # Record the start time of the command execution
        self._file_stat[job_id]["start timestamp"] = strftime("%Y-%m-%d %H:%M:%S", localtime())
        start_time = time()

        # Update the status of the job to "Running"
        self._file_stat[job_id]["status"] = "Running"
        self.to_csv(self._panel_name)

        # Execute the command
        try:
            self.run(cmd)
        finally:
            # Record the end time of the command execution
            self._file_stat[job_id]["end timestamp"] = strftime("%Y-%m-%d %H:%M:%S", localtime())
            end_time = time()

            # Calculate and record the duration of the command execution
            self._file_stat[job_id]["duration"] = round(end_time - start_time, 2)
    def run(self, cmd: str) -> None:
        """
        Execute a shell command.

        Args:
            cmd (str): The shell command to be executed.
        """
        import subprocess
        subprocess.run(cmd, shell=True)
    def parallel(self, processes: int = int(cpu_count() / 8)) -> None:
        """
        Refactored controller function to use ThreadPoolExecutor for parallel execution.
        
        Args:
            self: the instance of the class
            max_jobs (int): maximum number of jobs for ThreadPoolExecutor
        
        Returns:
            None
        """
        from concurrent.futures import ThreadPoolExecutor, as_completed
        self._parallel_jobs = processes
        
        # Create ThreadPoolExecutor with specified number of workers
        with ThreadPoolExecutor(processes) as executor:
            # Dictionary to store futures and corresponding job IDs
            futures = {}
            
            # Construct command for execution
            if self._cmds is None:
                self._cmds = []
                # Iterate over job dictionary
                for self._file, key in self._job_dict.items():
                    if self._format == "g16" or self._format == "ams" or self._format == "vasp":
                        self._default_cmd()
                    else:
                        raise ValueError("Unsupported format. Please use 'g16', or 'ams'.")
                    # Submit command to ThreadPoolExecutor and store future with job ID
                    futures[executor.submit(self.stamp_and_run, self._cmd, key)] = key
                    self._cmds.append(self._cmd)
            else:
                for idx, cmd in enumerate(self._cmds):
                    # Submit command to ThreadPoolExecutor and store future with job ID
                    futures[executor.submit(self.stamp_and_run, cmd, idx)] = idx
            
            # Iterate over completed futures
            for future in as_completed(futures):
                job_id = futures[future]
                try:
                    # Retrieve result from future
                    future.result()
                except Exception as e:
                    self._exception(job_id, e)
                else:
                    self._end(job_id)
                finally:
                    self.to_csv(self._panel_name)
    def build_panel(self, panel_name: str=".panel.csv") -> None:
        """
        Initialize the job queue and create a file status report.

        Args:
            format (str): format of the command

        Returns:
            None
        """
        self._move_to_queue(self._file_list)
        
        # Initialize a DataFrame to store file status
        self._file_stat = []

        # Iterate through the file list and populate the DataFrame and job dictionary
        for idx, file in enumerate(self._file_list):
            self._file_stat.append({
                "file": file.split('.')[0],
                "status": "None",
                "start timestamp": "None",
                "end timestamp": "None",
                "duration": "None",
                "error message": "None",
                "wrong message": "None"
            })
            self._job_dict[file.split('.')[0]] = idx
        self.to_csv(panel_name)
    def to_csv(self, panel_name: str=".panel.csv"):
        self._panel_name = panel_name
        # Save the file status report to a CSV file
        with open(self._panel_name, "w") as write_file:
            for idx, c in enumerate(self._file_stat[0].keys()):
                write_file.write(c if idx == 0 else f",{c}")
            write_file.write("\n")
            for r in self._file_stat:
                for idx, c in enumerate(r.values()):
                    write_file.write(c if idx == 0 else f",{c}")
                write_file.write("\n")
        with open("panel.csv", "w") as write_file:
            write_file.write(",".join(self._file_stat[0].keys()) + "\n")
            for panel_file in glob(".panel*.csv"):
                with open(panel_file, "r") as read_file:
                    for line in read_file.readlines()[1:]:
                        write_file.write(line)
    def _dir_check(self):
        """Check and create directories if they do not exist."""
        from os import makedirs
        for directory in ["done", "error", "wrong", "queue"]:
            makedirs(directory, exist_ok=True)

class Gaussian(Multiprocess):
    def __init__(self, file_list: list, **kwargs) -> None:
        super().__init__(file_list=file_list, **kwargs)
        self._format = "g16"
    def _default_cmd(self):
        self._cmd = f"cd {path.join(self._dir_path, 'queue')} && {self._format} < {self._file}.gjf > {self._file}.log"
    def _move_to_queue(self, file_list: list) -> None:
        super()._move_to_queue(file_list)
    def _exception(self, job_id: str, exception: str) -> None:
        super()._exception(job_id, exception)
        base_file = path.join(self._dir_path, "queue", self._file_stat[job_id]["file"])
        if path.isfile(base_file + ".log"):
            system(f"mv {base_file}.log {path.join(self._dir_path, 'error')}/")
    def _end(self, job_id: str) -> None:
        """
        Update the status of a job based on the content of its log file.

        This function will be called after the job is done. It will check the log file
        to see if the job has completed normally. If the log file contains
        "Normal termination", the status of the job will be updated to "Done"
        and the file will be moved to the "done" directory. Otherwise, the status
        of the job will be updated to "Wrong" and the file will be moved to the
        "wrong" directory.

        Args:
            job_id: str, the job ID in the first column of the panel.csv file

        Returns:
            None
        """

        # Get the file name based on the job_id
        base_file = path.join(self._dir_path, "queue", self._file_stat[job_id]["file"])
        log_path = base_file + ".log"
        queue_dir = path.join(self._dir_path, "queue")

        # Open the log file in read mode
        with popen(f"tail -10 {log_path}") as read_file:
            # Check if the log file contains "Normal termination"
            last_10_lines = read_file.readlines()[-10:]
            if "Normal termination" in last_10_lines[-1]:
                status = "done"
            else:
                status = "wrong"
                for idx, line in enumerate(last_10_lines):
                    if "termination" in line:
                        # Store the line before the line containing "termination" as the wrong message
                        self._file_stat[job_id]["wrong message"] = last_10_lines[idx-1].replace('\n', '')
                        break
            
            # Update the status to "Done" and move the file to the "done" directory
            self._file_stat[job_id]["status"] = status[0].upper() + status[1:]

            # Move the chk file to the "done" directory
            with open(f"{base_file}.gjf") as gjf_file:
                for line in gjf_file:
                    chk_file = line.split('=')[1].split()[0]
                    if "%chk" in line and path.isfile(path.join(queue_dir, chk_file)):
                        system(f"mv {path.join(queue_dir, chk_file)} {path.join(self._dir_path, status)}/")
                        break

            # Move the file to the "done" directory
            system(f"mv {base_file}.gjf {path.join(self._dir_path, status)}/")
            system(f"mv {base_file}.log {path.join(self._dir_path, status)}/")

class VASP(Multiprocess):
    def __init__(self, file_list: list, **kwargs) -> None:
        super().__init__(file_list=file_list, **kwargs)
        self._format = "vasp"
        self._partition = Partition()
    @property
    def processes(self):
        return self._partition.processes
    @processes.setter
    def processes(self, processes):
        self._partition.processes = processes
    def _default_cmd(self):
        self._cmd = f"cd {path.join(self._dir_path, 'queue', self._file)} && mpiexec.hydra -bootstrap fork -np {self._partition.processes} -ppn {self._partition.processes} vasp_std >& {self._file}.log"
    def _move_to_queue(self, file_list: list) -> None:
        from getpass import getuser
        for file in file_list:
            system(f"cp /home/{getuser()}/bin/vdw_kernel.bindat {path.join(self._dir_path, file)}")
        super()._move_to_queue(file_list)
    def _exception(self, job_id: str, exception: str) -> None:
        super()._exception(job_id, exception)
        system(f"mv {base_file} {path.join(self._dir_path, 'error')}/")
    def _end(self, job_id: str) -> None:
        """
        Update the status of a job based on the content of its log file.
        This function will be called after the job is done. It will check the log file
        to see if the job has completed normally. If the log file contains
        "General timing and accounting informations", the status of the job will be updated to "Done"
        and the file will be moved to the "done" directory. Otherwise, the status
        of the job will be updated to "Wrong" and the file will be moved to the
        "wrong" directory.
        Args:
            job_id: str, the job ID in the first column of the panel.csv file
        Returns:
            None
        """
        # Get the file name based on the job_id
        base_file = path.join(self._dir_path, "queue", self._file_stat[job_id]["file"])
        # Open the log file in read mode
        with popen(f"tail -20 {path.join(base_file, 'OUTCAR')}") as read_file:
            # Check if the log file contains "General timing and accounting informations"
            for line in read_file.readlines():
                if "General timing and accounting informations" in line:
                    # Update the status to "Done" and move the file to the "done" directory
                    self._file_stat[job_id]["status"] = "Done"
                    status = "done"
                    break
                elif any([word in line for word in ["ERROR", "internal error", "VERY BAD NEWS", "stop", "fatal"]]):
                    self._file_stat[job_id]["status"] = "Wrong"
                    self._file_stat[job_id]["wrong message"] = line.replace('\n', '')
                    status = "wrong"
                    break
            else:
                self._file_stat[job_id]["status"] = "Wrong"
                self._file_stat[job_id]["wrong message"] = "Unknown"
                status = "wrong"
        system(f"mv {base_file} {path.join(self._dir_path, status)}/")
