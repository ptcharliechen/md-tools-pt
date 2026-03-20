from kit.args import args
arg = args({"output": "rdf", "cpu": [-1, int, "specify the maximum number of applied CPU cores"], \
            "samename": [1, int, "0: without same name check, 1: with same name check"], "cn": [1, int, "0: no coordination number, 1: coordination number"], \
            "rmin": [0.0, float, "minimum distance"], "rmax": [10.0, float, "maximum distance"], "nbin": [501, int, "number of intervals"], \
            "plot": [None, int, "0: no plot, 1: combine RDF and CN, 2: do not combine RDF and CN"]}, "XDATCAR", molecules=True, output_isdir=True)
from os import path, mkdir, cpu_count
from re import compile
from concurrent.futures import ProcessPoolExecutor
from numpy import array_split, linspace
import numpy as np
from pandas import DataFrame
from kit.fundamental import AtomStep_Trj
from kit.function import RDF
from kit.vasp import XDATCAR
from kit.interface import blocks, step_lines

def function(num, total_AtomStep, center_AtomStep, measure_AtomStep, r_min, r_max, nbin):
    from scipy.interpolate import interp1d
    print(f"Calculating RDF {num}...")
    rdf = RDF()
    rdf.atom_step, rdf.center_atom_step, rdf.measure_atom_step = total_AtomStep, measure_AtomStep, center_AtomStep
    rdf_data, cn_data = rdf.run(r_min=r_min, r_max=r_max, nbins=nbin)
    rdf_fun = interp1d(linspace(r_min, r_max, nbin), rdf_data, kind="cubic")
    cn_fun = interp1d(linspace(r_min, r_max, nbin), cn_data, kind="cubic")
    rdf_df = DataFrame(columns=["interval", "rdf", "cn"])
    rdf_df["interval"] = np.round(linspace(r_min, r_max, nbin*10), 4)
    rdf_df["rdf"] = [rdf_fun(num)*rdf.center_position.shape[0] if rdf_fun(num) > 0 else 0 for num in rdf_df["interval"]]
    rdf_df["cn"] = [cn_fun(num)*rdf.center_position.shape[0] if cn_fun(num) > 0 else 0 for num in rdf_df["interval"]]
    return rdf_df

if __name__ == "__main__":
    xdatcar = XDATCAR(arg)

    atom_list, atom_list_flatten, atom_info = blocks(xdatcar, 2, flatten_flag=True, print_output_flag=True)
    step_list, step_list_flatten, step_info = step_lines(xdatcar, items=len(atom_list), flatten_flag=True, print_output_flag=True)

    if xdatcar.args.cpu <= 0:
        if cpu_count() < 4:
            xdatcar.args.cpt = 1
        elif len(step_list) > cpu_count()*0.25:
            xdatcar.args.cpu = int(cpu_count()*0.25)
        else:
            xdatcar.args.cpu = len(step_list)
    elif xdatcar.args.cpu > cpu_count():
        xdatcar.args.cpu = cpu_count()
    
    while(1):
        inp = input("Are coordination numbers needed? (y/n) [y]: ")
        if inp == '' or inp.lower() == 'y':
            cn_flag = True
            break
        elif inp.lower() == 'n':
            cn_flag = False
            break
        else:
            print("Warning: Input error.")
    while(1):
        inp = input("Minimum radius [0.0]: ")
        if inp == '':
            rmin = 0.0
            break
        elif compile(r"\d+\.?\d*").match(inp) is not None:
            rmin = float(inp)
            break
        else:
            print("Warning: Input error.")
    while(1):
        inp = input("Maximum radius [5.05]: ")
        if inp == '':
            rmax = 5.05
            break
        elif compile(r"\d+\.?\d*").match(inp) is not None:
            rmax = float(inp)
            break
        else:
            print("Warning: Input error.")
    while(1):
        inp = input("Number of slices [101]: ")
        if inp == '':
            nbin = 101
            break
        elif inp.isdigit():
            nbin = int(inp)
            if nbin % 5 == 0 or nbin % 10 == 0:
                rmax += (rmax-rmin) / nbin
                nbin += 1
            break
        else:
            print("Warning: Input error.")

    while(1):
        plot_flag, twin_flag = False, False
        inp = input("Plot (y/n) [n]: ")
        if inp.lower() == '':
            plot_flag = False
            break
        elif inp.lower() == 'y':
            plot_flag = True
            if cn_flag:
                while(1):
                    inp_2 = input("Combine rdf and coordination numbers in a figure (y/n) [y]: ")
                    if inp_2.lower() == '':
                        inp_2 = 'y'
                    if inp_2.lower() == 'y':
                        twin_flag = True
                        from kit.plot import TwinPlot
                        twin_plots = []
                        for i in range(len(atom_list)):
                            twin_plots.append(TwinPlot())
                        break
                    elif inp_2.lower() == 'n':
                        twin_flag = False
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
        elif inp.lower() == 'n':
            break
        else:
            print("Warning: Input error.")
    
    interval = np.round(linspace(rmin, rmax, nbin*10), 4)
    scripts, paras = [], []
    nums = []
    print("\nReading XDATCAR...")
    xdatcar.read_all(steps=step_list_flatten, atoms=atom_list_flatten)

    center_AtomSteps, measure_AtomSteps = [], []
    if xdatcar.args.cpu < len(step_list)+1:
        new_cores = [1 for _ in range(len(step_list))]
    else:
        weights = np.zeros(len(step_list))
        for idx, atoms in enumerate(atom_list):
            atoms_1, atoms_2 = atoms[0].get(), atoms[1].get()
            weights[idx] = len(atoms_1) * len(atoms_2)
        weights /= sum(weights)
        cores = weights * xdatcar.args.cpu
        new_weights = [weights[i] if cores[i] > 1 else 0 for i in range(len(weights))]
        new_weights /= sum(new_weights)
        remaining_cpu = xdatcar.args.cpu - sum([1 if cores[i] <= 1.0 else 0 for i in range(len(cores))])
        new_cores = [int(np.ceil(new_weights[i]*remaining_cpu)) if cores[i] > 1.0 else 1 for i in range(len(cores))]
        if sum(new_cores) > xdatcar.args.cpu:
            while(sum(new_cores) != xdatcar.args.cpu):
                new_cores[new_cores.index(max(new_cores))] -= 1
        elif sum(new_cores) < xdatcar.args.cpu:
            while(sum(new_cores) != xdatcar.args.cpu):
                new_cores[new_cores.index(min(new_cores))] += 1
        for i in range(len(new_cores)):
            if len(step_list[i].get()) < new_cores[i]:
                new_cores[i] = 1

    for idx, (atoms, steps) in enumerate(zip(atom_list, step_list)):
        atoms_1, atoms_2 = atoms[0].get(), atoms[1].get()
        step = array_split(steps.get(), new_cores[idx])
        for i in range(new_cores[idx]):
            print(f"Assigning steps {idx+1}_{i+1}...")
            center_AtomStep, measure_AtomStep = AtomStep_Trj(xdatcar), AtomStep_Trj(xdatcar)
            center_AtomStep.lattice = measure_AtomStep.lattice = xdatcar.lattice
            center_AtomStep.atoms, measure_AtomStep.atoms = atoms[0], atoms[1]
            center_AtomStep.steps.put(step[i])
            center_AtomStep.split(xdatcar.atom_step, center_AtomStep.steps, atoms[0])
            measure_AtomStep.split(xdatcar.atom_step, center_AtomStep.steps, atoms[1])
            center_AtomSteps.append(center_AtomStep); measure_AtomSteps.append(measure_AtomStep)
            nums.append(f"{idx+1}_{i+1}")
    
    try:
        with ProcessPoolExecutor(max_workers=sum(new_cores)) as executor:
            futures = [executor.submit(function, num, xdatcar.atom_step, center_AtomStep, measure_AtomStep, rmin, rmax, nbin) for num, center_AtomStep, measure_AtomStep in zip(nums, center_AtomSteps, measure_AtomSteps)]
            results = [future.result() for future in futures]
    except Exception as err:
        print(err)
        print("Error: Some troubles happened during parallel calculation. Calculate again.")
        exit()
    line_num = idx + 1
    
    rdf_array = np.zeros((len(step_list)+1, len(interval)))
    rdf_array[0] = interval
    if cn_flag:
        cn_array = np.zeros((len(step_list)+1, len(interval)))
        cn_array[0] = interval
    
    for num, result in zip(nums, results):
        idx = int(num.split("_")[0])
        rdf_array[idx] += result["rdf"]
        if cn_flag:
            cn_array[idx] += result["cn"]
    for i in range(line_num):
        rdf_array[i+1] /= len(step_list[i].get())
        if cn_flag:
            cn_array[i+1] /= len(step_list[i].get())
    
    if twin_flag or plot_flag:
        for num in nums:
            idx = int(num.split("_")[0])
            if twin_flag:
                twin_plots[idx-1].append(interval, list(rdf_array[idx]), axis="left")
                twin_plots[idx-1].append(interval, list(cn_array[idx]), axis="right")
                twin_plots[idx-1].xlabel = "Radius (Å)"
                twin_plots[idx-1].left_ylabel = "g(r)"
                twin_plots[idx-1].right_ylabel = "CN"
                twin_plots[idx-1].left_ylim = (-0.01, None)
                twin_plots[idx-1].right_ylim = (-0.001, None)
            elif plot_flag:
                rdf_plot.append(interval, list(rdf_array[idx]))
                rdf_plot.xlabel = "Radius (Å)"
                rdf_plot.ylabel = "g(r)"
                rdf_plot.ylim = (-0.01, None)
                if cn_flag:
                    cn_plot.append(interval, list(cn_array[idx]))
                    cn_plot.xlabel = "Radius (Å)"
                    cn_plot.ylabel = "CN"
                    cn_plot.ylim = (-0.001, None)
    
    rdf_df = DataFrame(rdf_array.T, columns=["interval"]+list(range(1, line_num+1)))
    cn_df = DataFrame(cn_array.T, columns=["interval"]+list(range(1, line_num+1)))

    atom_indices = "atoms,"+"".join([f",{atom_info[i-1][0].replace(',', ';')} {atom_info[i-1][1].replace(',', ';')}" for i in range(1, line_num+1)])
    step_indices = "steps,"
    for i in range(line_num):
        step_indices += ","
        if str(step_info[i]).lower() != "all":
            for idx, step in enumerate(step_info[i]):
                if idx > 0:
                    step_indices += f"_{step.replace(',', ';')}"
                else:
                    step_indices += f"{step.replace(',', ';')}"
        else:
            step_indices += step_info[i].lower()
    mkdir(arg.args.output)
    rdf_df.to_csv(path.join(arg.args.output, "rdf.csv"), float_format="%.4f")
    with open(path.join(arg.args.output, "rdf.csv")) as read_file:
        lines = read_file.readlines()
        lines.insert(1, atom_indices+"\n")
        lines.insert(2, step_indices+"\n")
    with open(path.join(arg.args.output, "rdf.csv"), "w") as write_file:
        for line in lines:
            write_file.write(line)
    if cn_flag:
        cn_df.to_csv(path.join(arg.args.output, "cn.csv"), float_format="%.4f")
        with open(path.join(arg.args.output, "cn.csv")) as read_file:
            lines = read_file.readlines()
            lines.insert(1, atom_indices+"\n")
            lines.insert(2, step_indices+"\n")
        with open(path.join(arg.args.output, "cn.csv"), "w") as write_file:
            for line in lines:
                write_file.write(line)
    
    if twin_flag or plot_flag:
        from kit.plot import plot
        if len(step_list) > 1:
            palette = None
            while(1):
                inp = input("Unify the palette settings (y/n) [y]: ")
                if inp.lower() == "":
                    inp = "y"
                if inp.lower() == "n":
                    break
                elif inp.lower() == "y":
                    break
                else:
                    print("Warning: Input error.")
        elif len(step_list) == 1:
            inp = "y"
    
        if inp.lower() == "y":
            while(1):
                inp_2 = input("Adjust color (y/n) [n]: ")
                if inp_2 == '' or inp_2.lower() == 'n':
                    from seaborn import color_palette
                    palette = color_palette("Set2", 8)
                    if len(step_list) > 8:
                        palette += color_palette("husl", len(step_list)-8)
                    break
                palette = []
                for i in range(len(step_list)):
                    inp_3 = input(f"Color {i+1}: ")
                    palette.append(inp_3)
                if twin_flag:
                    for idx in range(len(atom_list)):
                        twin_plots[idx].palette = palette
                elif plot_flag:
                    rdf_plot.palette = palette
                    if cn_flag:
                        cn_plot.palette = palette
                break

    if plot_flag:
        mkdir(path.join(arg.args.output, "fig"))
    if twin_flag:
        for idx in range(1, len(atom_list)+1):
            if palette is None:
                print(f"Plot {idx}")
            plot(twin_plots[idx-1], path.join(arg.args.output, "fig", f"{idx}.png"), line_label=False, palette=palette)
    elif not twin_flag and plot_flag:
        print("RDF Plot")
        plot(rdf_plot, path.join(arg.args.output, "fig", "rdf.png"), palette=palette)
        if cn_flag:
            print("CN Plot")
            if rdf_plot.line_label != False:
                cn_plot.line_label = rdf_plot.line_label
            plot(cn_plot, path.join(arg.args.output, "fig", "cn.png"), palette=palette)

    print("Done!")
