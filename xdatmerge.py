from os import path, getcwd
from sys import argv
from kit.fundamental import Args

if __name__ == "__main__":
    Args.same_name(getcwd(), "XDATCAR")

    i = 1
    XDAT_paths = []
    for idx, XDAT_path in enumerate(argv):
        if idx != 0:
            if path.isfile(XDAT_path):
                XDAT_paths.append(XDAT_path)
            elif path.isfile(path.join(XDAT_path, "XDATCAR")):
                XDAT_paths.append(path.join(XDAT_path, "XDATCAR"))
            else:
                print("Warning: Argument {} doesn't exist.".format(idx))
    if len(argv) == 1: 
        print("Input paths of XDATCAR.\n")
        while(1):
            XDAT_path = input("%d: " % (i))
            if XDAT_path.lower() == "end":
                break
            elif path.isfile(XDAT_path):
                i += 1
                XDAT_paths.append(XDAT_path)
            elif path.isfile(path.join(XDAT_path, "XDATCAR")):
                i += 1
                XDAT_paths.append(path.join(XDAT_path, "XDATCAR"))
            else:
                print("Warning: This file doesn't exist.")

    num = 0
    with open(path.join(getcwd(), "XDATCAR"), "w") as write_file:
        for idx, XDAT_path in enumerate(XDAT_paths):
            with open(XDAT_path) as read_file:
                row = 0
                for line in read_file:
                    if "Direct configuration" in line:
                        num += 1
                        if int(line.split()[-1]) != int(num):
                            line = line.replace(line[-6:-1], ("%6d" % (num)))
                    row += 1
                    if idx > 0 and row < 8:
                        pass
                    else:
                        write_file.write(line)
            if '\n' not in line:
                write_file.write('\n')
    print("Done!")
