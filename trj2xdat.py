from sys import argv
from os import path, getcwd, listdir
from kit.fundamental import Args
from kit.software import Conquest, QE
from kit.vasp import XDATCAR

if __name__ == "__main__":
    arg = Args({"input": getcwd(), "output": "XDATCAR", "elementlist": None})
    arg.same_name(getcwd(), "XDATCAR")
    if len(argv) == 1 or (len(argv) == 2 and not path.isfile(argv[1]) and not path.isdir(argv[1])):
        for f in listdir(getcwd()):
            if "mdtrj" in f:
                arg.args.input = path.join(getcwd(), f)
                arg.args.elementlist = path.join(getcwd(), "QE.in")
            elif "xsf" in f:
                arg.args.input = path.join(getcwd(), f)
    elif path.isfile(argv[1]):
        arg.args.input = argv[1]
    elif not path.isfile(argv[1]) and path.isdir(argv[1]):
        for f in listdir(argv[1]):
            if "mdtrj" in f:
                arg.args.input = path.join(argv[1], f)
            elif "xsf" in f:
                arg.args.input = path.join(argv[1], f)
    else:
        arg.input_file_check(getcwd(), argv[1])
    
    while(1):
        if "mdtrj" in arg.args.input:
            kind = 'q'
            break
        elif "xsf" in arg.args.input:
            kind = 'c'
            break
        else:
            kind = input("Software, QE or Conquest (q/c): ")
            if kind.lower() == 'q':
                arg.input_type = "QE"
            elif kind.lower() == 'c':
                arg.input_type = "Conquest"
            else:
                print("Warning: Input 'q' or 'c'.")
                continue
            kind = kind.lower()
            arg.input_check()
            break
    if kind == 'q':
        if path.isfile(arg.args.elementlist):
            pass
        elif path.isfile(path.join(arg.args.input, "QE.in")):
            arg.args.elementlist = path.join(arg.args.input, "QE.in")
        else:
            arg.args.elementlist = arg.input_file_check(arg.args.elementlist, "QE.in")
        qe = QE(arg)
        qe.read_all()
        xdatcar = XDATCAR(arg)
        xdatcar.bridge = qe
        xdatcar.write_all()
    elif kind == 'c':
        conquest = Conquest(arg)
        xdatcar = XDATCAR(arg)
        conquest.read_all()
        xdatcar.bridge = conquest
        xdatcar.write_all()
    print("Done!")
