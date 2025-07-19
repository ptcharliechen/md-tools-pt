from os import getcwd, cpu_count
from concurrent.futures import ProcessPoolExecutor
from numpy import array, array_split, linspace
import numpy as np
from pandas import DataFrame
from scipy.interpolate import interp1d
from kit.fundamental import AtomStep_Trj
from kit.function import RDF
from kit.vasp import XDATCAR
from kit.interface import args, blocks, step_lines
from kit.plot import plot
from re import compile

def function(num, total_AtomStep, center_AtomStep, measure_AtomStep, r_min, r_max, nbin):
    print(f"Calculating RDF {num}...")
    rdf = RDF()
    rdf.atom_step, rdf.center_atom_step, rdf.measure_atom_step = total_AtomStep, measure_AtomStep, center_AtomStep
    rdf.run(r_min=r_min, r_max=r_max, nbins=nbin)
    #rdf_fun = interp1d(linspace(r_min, r_max, nbin), rdf.rdf, kind="cubic")
    #cn_fun = interp1d(linspace(r_min, r_max, nbin), rdf.coordination_number, kind="cubic")
    rdf_df = DataFrame(columns=["interval", "rdf", "cn"])
    rdf_df["interval"] = rdf.interval[:-1]
    rdf_df["rdf"] = rdf.rdf * len(rdf.center_position)
    rdf_df["cn"] = rdf.coordination_number * len(rdf.center_position)
    return rdf_df

if __name__ == "__main__":
    arg = args({"output": "rdf", "molecules": [None, str, "specify a file with user-interested molecules"], \
                "cpu": [int(cpu_count()/2.5), int, "specify the maximum number of applied CPU weights"]}, "XDATCAR")
    arg.same_name(getcwd(), "atom_list.csv")
    
    xdatcar = XDATCAR(arg)
    
    atom_list, atom_list_flatten = blocks(xdatcar, 2, flatten_flag=True)
    step_list, step_list_flatten = step_lines(xdatcar, items=len(atom_list), flatten_flag=True)
    step_list_idxs = []
    for step_list_single in step_list:
        step_list_idx = []
        for step in step_list_single.get():
            step_list_idx.append(step_list_flatten.index_list[step])
        step_list_idxs.append(step_list_idx)
    while(1):
        tmp = input("Are coordination numbers needed? (y/n) [y]: ")
        if tmp == '' or tmp.lower() == 'y':
            cn_flag = True
            arg.same_name(getcwd(), "cn.csv")
            break
        elif tmp.lower() == 'n':
            cn_flag = False
            break
        else:
            print("Warning: Input error.")
    while(1):
        tmp = input("Minimum radius [0.0]: ")
        if tmp == '':
            rmin = 0.0
            break
        elif compile(r"\d+\.?\d*").match(tmp) is not None:
            rmin = float(tmp)
            break
        else:
            print("Warning: Input error.")
    while(1):
        tmp = input("Maximum radius [5.05]: ")
        if tmp == '':
            rmax = 5.05
            break
        elif compile(r"\d+\.?\d*").match(tmp) is not None:
            rmax = float(tmp)
            break
        else:
            print("Warning: Input error.")
    while(1):
        tmp = input("Number of slices [101]: ")
        if tmp == '':
            nbin = 101
            break
        elif tmp.isdigit():
            nbin = int(tmp)
            if nbin % 5 == 0 or nbin % 10 == 0:
                rmax += (rmax-rmin) / nbin
                nbin += 1
            break
        else:
            print("Warning: Input error.")

    while(1):
        plotFlag, twinFlag = False, False
        tmp = input("Plot (y/n) [n]: ")
        if tmp.lower() == '':
            plotFlag = False
            break
        elif tmp.lower() == 'y':
            arg.same_name(getcwd(), f"{arg.args.output}.png")
            plotFlag = True
            if cn_flag:
                while(1):
                    tmp = input("Combine rdf and coordination numbers in a figure (y/n) [y]: ")
                    if tmp.lower() == '':
                        twinFlag = True
                        from kit.plot import TwinPlot
                        twin_plots = []
                        for i in range(len(atom_list)):
                            twin_plots.append(TwinPlot())
                        break
                    elif tmp.lower() == 'y':
                        twinFlag = True
                        from kit.plot import TwinPlot
                        twin_plots = []
                        for i in range(len(atom_list)):
                            twin_plots.append(TwinPlot())
                        break
                    elif tmp.lower() == 'n':
                        twinFlag = False
                        from kit.plot import Plot
                        rdf_plot = Plot()
                        cn_plot = Plot()
                        break
                    else:
                        print("Warning: Input error.")
            else:
                from kit.plot import Plot
                rdf_plot = Plot()
            break
        elif tmp.lower() == 'n':
            break
        else:
            print("Warning: Input error.")
    
    #interval = np.round(linspace(rmin, rmax, nbin*10), 4)
    scripts, paras = [], []
    nums = []
    print("\nReading XDATCAR...")
    xdatcar.read_all(steps=step_list_flatten, atoms=atom_list_flatten)

    center_AtomSteps, measure_AtomSteps = [], []
    if arg.args.cpu < len(step_list)+1:
        NewCores = [1 for _ in range(len(step_list))]
    else:
        weights = np.zeros(len(step_list))
        for idx, atoms in enumerate(atom_list):
            atoms_1, atoms_2 = atoms[0].get(), atoms[1].get()
            weights[idx] = len(atoms_1) * len(atoms_2)
        weights /= sum(weights)
        cores = weights * arg.args.cpu
        NewWeights = [weights[i] if cores[i] > 1 else 0 for i in range(len(weights))]
        NewWeights /= sum(NewWeights)
        RemainingCPU = arg.args.cpu - sum([1 if cores[i] <= 1.0 else 0 for i in range(len(cores))])
        NewCores = [int(np.ceil(NewWeights[i]*RemainingCPU)) if cores[i] > 1.0 else 1 for i in range(len(cores))]
        if sum(NewCores) > xdatcar.args.cpu:
            while(sum(NewCores) != xdatcar.args.cpu):
                NewCores[NewCores.index(max(NewCores))] -= 1
        elif sum(NewCores) < xdatcar.args.cpu:
            while(sum(NewCores) != xdatcar.args.cpu):
                NewCores[NewCores.index(min(NewCores))] += 1
        for i in range(len(NewCores)):
            if len(step_list[i].get()) < NewCores[i]:
                NewCores[i] = 1

    for idx, (atoms, steps) in enumerate(zip(atom_list, step_list)):
        atoms_1, atoms_2 = atoms[0].get(), atoms[1].get()
        step = array_split(steps.get(), NewCores[idx])
        cen_pos, mea_pos = [], []
        for i, atom_idx in enumerate(atom_list_flatten.get()):
            if atom_idx in atoms_1:
                cen_pos.append(xdatcar.fractional_position[step_list_idxs[idx], i])
            elif atom_idx in atoms_2:
                mea_pos.append(xdatcar.fractional_position[step_list_idxs[idx], i])
        cen_pos = array_split(array(cen_pos).transpose((1, 0, 2)), NewCores[idx])
        mea_pos = array_split(array(mea_pos).transpose((1, 0, 2)), NewCores[idx])
        for i in range(NewCores[idx]):
            print(f"Assigning steps {idx+1}_{i+1}...")
            center_AtomStep, measure_AtomStep = AtomStep_Trj(xdatcar), AtomStep_Trj(xdatcar)
            center_AtomStep.lattice = measure_AtomStep.lattice = xdatcar.lattice
            center_AtomStep.atoms, measure_AtomStep.atoms = atoms[0], atoms[1]
            center_AtomStep.steps.put(step[i]); measure_AtomStep.steps.put(step[i])
            center_AtomStep.fractional_position, measure_AtomStep.fractional_position = cen_pos[i], mea_pos[i]
            center_AtomSteps.append(center_AtomStep); measure_AtomSteps.append(measure_AtomStep)
            nums.append(f"{idx+1}_{i+1}")
    
    try:
        with ProcessPoolExecutor(sum(NewCores)) as executor:
            futures = [executor.submit(function, num, xdatcar.atom_step, center_AtomStep, measure_AtomStep, rmin, rmax, nbin) for num, center_AtomStep, measure_AtomStep in zip(nums, center_AtomSteps, measure_AtomSteps)]
            results = [future.result() for future in futures]
    except:
        print("Error: Some troubles happened during parallel calculation. Calculate again.")
        exit()
    
    LineNum = idx + 1
    interval = results[0]["interval"]
    rdf_array = np.zeros((len(step_list)+1, len(interval)))
    rdf_array[0] = interval
    atom_list_df = DataFrame(columns=["i", "j"])
    atom_list_df.index += 1
    if cn_flag:
        cn_array = np.zeros((len(step_list)+1, len(interval)))
        cn_array[0] = interval
    print(cn_array.shape)
    for i in range(LineNum):
        atom_list_df.loc[i+1, "i"] = str(atom_list[i][0].get()).replace('(', '').replace(')', '').replace(',', '')
        atom_list_df.loc[i+1, "j"] = str(atom_list[i][1].get()).replace('(', '').replace(')', '').replace(',', '')
    for num, result in zip(nums, results):
        idx = int(num.split("_")[0])
        rdf_array[idx] += result["rdf"]
        if cn_flag:
            cn_array[idx] += result["cn"]
        if twinFlag:
            twin_plots[idx-1].append(interval, list(rdf_array[idx]), axis="left")
            twin_plots[idx-1].append(interval, list(cn_array[idx]), axis="right")
            twin_plots[idx-1].xlabel = "Radius (Å)"
            twin_plots[idx-1].left_ylabel = "g(r)"
            twin_plots[idx-1].right_ylabel = "CN"
        elif plotFlag:
            rdf_plot.append(interval, list(rdf_array[idx]))
            rdf_plot.xlabel = "Radius (Å)"
            rdf_plot.ylabel = "g(r)"
            if cn_flag:
                cn_plot.append(interval, list(cn_array[idx]))
                cn_plot.xlabel = "Radius (Å)"
                cn_plot.ylabel = "CN"
    interval = np.round(linspace(rmin, rmax, nbin*10), 4)
    rdf_df_array = np.zeros((len(step_list)+1, len(interval)))
    if cn_flag:
        cn_df_array = np.zeros((len(step_list)+1, len(interval)))
    for idx in range(LineNum):
    rdf_df_array
    rdf_df = DataFrame(rdf_array.T, columns=["interval"]+list(range(1, LineNum+1)))
    cn_df = DataFrame(cn_array.T, columns=["interval"]+list(range(1, LineNum+1)))
    for i in range(LineNum):
        rdf_df[i+1] /= len(step_list[i].get())
        if cn_flag:
            cn_df[i+1] /= len(step_list[i].get())
    
    atom_list_df.to_csv("atom_list.csv")
    rdf_df.to_csv(arg.args.output+".csv", float_format="%.4f")
    if cn_flag:
        cn_df.to_csv("cn.csv", float_format="%.4f")
    if twinFlag:
        for idx in range(1, len(atom_list)+1):
            print(f"Plot {idx}")
            arg.same_name(getcwd(), f"{idx}.png")
            plot(twin_plots[idx-1], f"{idx}.png", line_label=False)
    elif plotFlag:
        print("RDF Plot")
        plot(rdf_plot, arg.args.output+".png")
        if cn_flag:
            print("CN Plot")
            plot(cn_plot, "cn.png")
    print("Done!")
