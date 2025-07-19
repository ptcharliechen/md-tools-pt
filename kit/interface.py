from os import path
from re import compile
from collections.abc import Iterable
from copy import deepcopy
from numpy import inf
from kit.fundamental import Step, Atom, Args

def args(parameters: dict, inputType: str, outputFile=None) -> Args:
    """
    Check and set the input and output files for preprocessing.

    input check:
        The input file exists or not. If not, reset the path.
    output check:
        The output file exists or not. If so, revise the old file.

    Parameters
    ----------
    parameters (mandatory): dict
        A dictionary of the parameter options allowed by the designer.
    inputType (mandatory): str
        Filename extension or default input file name. Choose one.
        'input_file_check' method helps you judge it.
        Attention, default input file name has higher priority.
        If not, filename extension is applied.
        Default input file name:
            If users just provide a directory, 'input_file_check' method implements it
            with 'inputType'.
        Filename extension:
            If users give a file without filename extension, 'input_file_check' method
            supplements it with 'inputType'.
    outputFile (optional): str, default None
        The default output file name. If None, the output file name will be the same as
        the 'output' key in 'parameters' dictionary.

    Returns
    -------
    args : Args
        A dictionary of preprocessing options.
    """
    # Initialize Args with parameters
    args = Args(parameters)
    
    # Check the 'input' key in 'parameters' dictionary is complete or not.
    # If not, a valid file path should be reset to 'input' key in 'args' by users.
    args.args.input = args.input_file_check(args.args.input, inputType)
    
    # Check the 'output' key in 'parameters' dictionary exists any same name file or directory.
    # If so, the old file or directory name should be revised by users.
    from os import getcwd
    if outputFile is None:
        methods = [attr for attr in dir(args.args) if "output" in attr]
        for method in methods:
            if getattr(args.args, method) is None:
                continue
            args.same_name(getcwd(), getattr(args.args, method))
    else:
        args.same_name(getcwd(), outputFile)
    
    return args

def flatten(item_list: list):
    """
    Flatten a list of Atom or Step objects into a single list of Atom objects.

    Parameters
    ----------
        item_list (mandatory): Iterable
            A list of Atom or Step objects.

    Returns
    -------
        Atom or Step :
            An Atom or Step object with the flattened list of Atom objects assigned to it.
    
    """
    def recursive(item_list):
        flatten = []
        for item in item_list:
            if isinstance(item, Iterable):
                _, item_object = recursive(item)
                flatten.extend(_)
            elif hasattr(item, 'get'):
                if isinstance(item, Atom):
                    item_list = list(item.get())
                    item_object = Atom(item.args)
                else:
                    item_list = list(item.get(slice_flatten=True))
                    item_object = Step(item)
                flatten.extend(item_list)
            else:
                raise ValueError("'item_list' should be a list of 'Atom', 'Step', or lists of 'Atom' or 'Step' objects.")
        return flatten, item_object
    
    flatten_list, item_object = recursive(item_list)
    flatten_list = list(set(flatten_list))
    
    item_object.put(flatten_list)

    # Sort the Atom objects in item_object if item_object is an Atom object
    if isinstance(item_object, Atom):
        item_object.sort()
    
    return item_object

def step(software: object, items: int=-1, information: bool=True, step=Step(), annotate=0) -> Step:
    """
    Prompt the user to input step numbers. All items will be combined into a 'Step' object.

    The user will be asked to input the step numbers one by one, until "end"
    is input.
    
    If "all" is input, all steps from the start step to the last step will be selected.
    
    Negative numbers mean counting from the last.
    
    The input can also be a range of step numbers, separated by a colon ":".

    Parameters
    ----------
        software (mandatory): object, no default
            The software object with the necessary methods and attributes.
        
        information (optional) : bool, default is True
            If True, print the information of start and step number.
        
        items (required): int, default is -1
            The number of steps to be input. The default is -1, which means input
            continues until the user input "end".

    Returns
    -------
        Step :
            A Step object that contains the step numbers.

    Examples
    --------
        "1" will select the first step.
        "-5" will select the last fifth step.
        "1:5" will select steps 1 to 5.
        "-5:-1" will select the last five steps.
    """
    # Initialize the step object and set its preprocess attribute
    i = 1
    step.bridge = software
    
    # Read step numbers from a file
    step_list_from_file = file(software.args.stepfile, items)
    
    # Print the start step number
    print(f"\nStart step number: {software.args.step}")
    if step_list_from_file is None and information:
        print("Input the numbers of steps.\n")

    # Loop until the user input "end" or the specified number of steps is reached
    while(1):
        if items > 0 and i > items:
            break
        
        if step_list_from_file is None or step_list_from_file == []:
            inp = input("%d: " % (i))
        else:
            inp = step_list_from_file.pop(0)
            print(f"{i}: {inp}")
        
        if annotate:
            inp, annotate_info = inp.split()[0], inp.split()[1:1+annotate]

        # Put the user input into the step object and check for errors
        err_return = step.put(inp)
        if step.err_list[err_return]:
            if step_list_from_file is not None:
                print(f"Line {i} in {software.args.stepfile}: {step.err_list[err_return]}")
            else:
                print(step.err_list[err_return])
            continue

        # Break the loop if "end" is input and the number of steps is specified
        if err_return == 1:
            if items > -1:
                print("Warning: 'end' is invalid when the line number is specified.")
                continue
            else:
                break
        # Break the loop if a range of steps is input
        elif err_return == 2:
            return "all"

        # Increment the step counter
        i += 1

    # Return the step object
    if annotate:
        return step, annotate_info
    else:
        return step

def step_lines(software: object, step=Step(), items: int=-1, information: bool=True, flatten_flag: bool=False, annotate: int=0) -> tuple:
    """
    Prompt the user to input step numbers. Every item will be a 'Step' object individually.
    
    The user will be asked to input the step numbers, until "end" is input.
    
    If "all" is input, all steps from the start step to the last step will be selected.
    
    If "same" is input, copy the last item.
    
    If users would like to input step numbers without continuation, separated by a camma ',' (See the last example).
    
    Negative numbers mean counting from the last.
    
    The input can also be a range of step numbers, separated by a colon ":".
    
    If flatten_flag is True, return an additional 'Step' object which contains all non-repetitive step numbers which users involve in.
    
    Parameters
    ----------
        software (mandatory) : object, no default
            The software object with the necessary methods and attributes.

        items (required) : int, default is -1
            The number of steps to be input.

        information (optional) : bool, default is True
            If True, print the information of start and step number.

        flatten_flag (optional) : bool, default is False
            If True, return an additional and non-repetitive step numbers in a 'Step' object.
    
    Returns
    -------
        list if flatten_flag is False:
            Contain a list of 'Step' objects. Every object contain the step numbers every item.

        tuple containing one list and one 'Step' object if flatten_flag is True:
            The first parameter as the previous list, and the second one contains the flatten list of all steps.
            It is a useful parameter if users would like to use 'read_step' method.
    
    Examples
    --------
        "1" will select the first step.
        "-5" will select the last fifth step.
        "1:5" will select steps 1 to 5.
        "-5:-1" will select the last five steps.
        "1:5,11:15" will select steps 1 to 5 and 11 to 15.
    """
    
    # Initialize
    i = 1
    all_flag = False
    step_list = []  # List of step numbers
    
    # Read step numbers from a file
    step_list_from_file = file(software.args.stepfile, items)
    
    # Print the start step number
    print(f"\nStart step number: {software.args.step}")
    if step_list_from_file is None and information:
        print("Input the numbers of steps.\n")
    
    pre_steps = deepcopy(step)
    pre_steps.bridge = software
    
    # Prompt the user to input step numbers
    annotate_info = []
    while(1):
        if items > 0 and i > items:
            # Add steps without repetation to the flattened list if "flatten" is True
            if flatten_flag:
                if annotate:
                    return (step_list, flatten(step_list), annotate_info) if not all_flag else (step_list, "all", annotate_info)
                else:
                    return (step_list, flatten(step_list)) if not all_flag else (step_list, "all")
            else:
                if annotate:
                    return step_list, annotate_info
                else:
                    return step_list
        
        if step_list_from_file is None or step_list_from_file == []:
            inp = input("%d: " % (i))
        else:
            inp = step_list_from_file.pop(0)
            print(f"{i}: {inp}")
        
        # Create a Step object
        steps = deepcopy(step)
        steps.bridge = pre_steps
        
        # Check if the input is "end"
        if inp.lower() == "end":
            if items > -1:
                print("Warning: 'end' is invalid when the number of the atom lines is specified.")
                continue
            elif flatten_flag:
                if annotate:
                    return step_list, flatten(step_list), annotate_info
                else:
                    return step_list, flatten(step_list)
            else:
                if annotate:
                    return step_list, annotate_info
                else:
                    return step_list
        # Check if the input is "all"
        elif inp.split()[0].lower() == "all":
            _ = steps.put(inp)
            all_flag = True
            step_list.append(steps)
            i += 1
            continue
        # Check if the input is "same"
        elif inp.split()[0].lower() == "same":
            step_list.append(step_list[-1])
            i += 1
            continue
        elif ":" in inp.lower() and "," in inp.lower():
            print("Warning: ':' and ',' can not be used at the same time.")
            continue
        
        if annotate:
            inp, tmp = inp.split()[0], inp.split()[1:1+annotate]
        sp = inp.split(',')
        
        # Check each step range input
        for step_range in sp:
            err_return = steps.put(step_range)
            if steps.err_list[err_return]:
                if step_list_from_file is not None:
                    print(f"Line {i} in {software.args.stepfile}: {steps.err_list[err_return]}")
                else:
                    print(steps.err_list[err_return])
                break
        else:
            i += 1
            step_list.append(steps)
            if annotate:
                annotate_info.append(tmp)

def single(software: object, atom_number: int=1, atom=Atom(), name: str="Atom", annotate: int=0) -> Atom:
    """
    Prompt user to input a single atom or a list of atoms.
    
    Args:
    -----
        software (mandatory): object, no default
            The software object with the necessary methods and attributes.
        atom_number (required): int, default is 1
            The number of atoms to expect as input.
        name (required): str, default is "Atom"
            The name to be displayed when requesting user input.
    
    Returns:
    --------
        list: An 'Atom' object representing the input atoms.
    """
    # Read the elements from the software object
    software.read_elements()
    
    # Read atom numbers from a file
    atom_list_from_file = file(software.args.atomfile, 1)
    
    # Display the start atom number if information is True
    print(f"\nStart atom number: {software.args.atom}\n")
    
    # Prompt the user to input atom numbers until a valid input is given
    while(1):
        # Prompt the user to input a single atom or a list of atoms
        if atom_list_from_file is None or atom_list_from_file == []:
            tmp = input(f"{name}: ")
        else:
            tmp = atom_list_from_file.pop(0)
            print(f"{name}: {tmp}")
        
        # Check the number of atoms input
        if atom_number == 1 and ' ' in tmp:
            print("Warning: One number is allowed only.")
            continue
        elif atom_number > 1 and ' ' not in tmp:
            print(f"Warning: {atom_number} number separated by {atom_number-1} space(s).")
            continue
        elif tmp == "":
            print("Warning: 'Null' is not allowed.")
            continue
        
        # Create an Atom object and set its bridge to the software object
        atoms = deepcopy(atom)
        atoms.bridge = software
        
        # Loop over the individual atoms in the input
        for a in tmp.split()[:atom_number]:
            # Check if the input is an integer
            if compile(r"^\-?\d+$").match(a) is None:
                print("Warning: Input integer(s).")
                break
            else:
                err_return = atoms.put(a)
                # Check if there is any error
                if atoms.err_list[err_return]:
                    print(atoms.err_list[err_return])
        else:
            if annotate:
                if atom_number < len(tmp.split()) <= atom_number + annotate:
                    return atoms, tmp.split()[atom_number:]
                elif len(tmp.split()) > atom_number + annotate:
                    return atoms, tmp.split()[atom_number:atom_number+annotate]
                else:
                    return atoms, []
            else:
                return atoms

def couple(software: object, block_number: int=2, items: int=-1, atom=Atom(), information: bool=True, Atom_split_flag: bool=False, flatten_flag: bool=False, annotate: int=0) -> list:
    """
    Prompt user to input a list of atoms.
    
    'couple' function only allows an atom every block, and 'blocks' function allows over one atom every block.
    
    Args:
    -----
        software (mandatory): object, no default
            The software object with the necessary methods and attributes.
        block_number (required): int, default is 2
            The number of atoms to expect as input per step.
        items (required): int, default is -1
            The number of steps to expect as input. If -1, the input will continue until "end" is input.
        information (optional): bool, default is True
            Whether to display information about the start atom number.
    
    Returns:
    --------
        list: A list of 'Atom' objects.
            Every 'Atom' object contains atom number whose number is equal to 'block_number'.
            Elements in an 'Atom' object are not sorted.
    """
    
    # Initialize counter
    i = 1
    
    # Read the elements from the software object
    software.read_elements()
    
    # Initialize list to store the input atoms
    atom_list = []
    
    # Read atom numbers from a file
    atom_list_from_file = file(software.args.atomfile, items)
    
    # Display the start atom number if information is True
    print(f"\nStart atom number: {software.args.atom}")
    if atom_list_from_file is None and information:
        print(f"Input {block_number} numbers of atoms separated by {block_number-1} space(s).\n")
    
    if annotate:
        annotate_info = []

    # Prompt the user to input atom numbers
    # Loop until the user input "end" or the specified number of steps is reached
    while(1):
        # Check if the number of steps is reached
        if items > 0 and i > items:
            if flatten_flag:
                if annotate:
                    return atom_list, flatten(atom_list), annotate_info
                else:
                    return atom_list, flatten(atom_list)
            else:
                if annotate:
                    return atom_list, annotate_info
                else:
                    return atom_list
        
        # Prompt the user to input a list of atoms
        if atom_list_from_file is None or atom_list_from_file == []:
            inp = input("%d: " % (i))
        else:
            inp = atom_list_from_file.pop(0)
            print(f"{i}: {inp}")
        
        # Check if the input is "end"
        if inp.lower() == "end":
            if items > -1:
                print("Warning: 'end' is invalid when the number of the step lines is specified.")
                continue
            else:
                if flatten_flag:
                    if not annotate:
                        return atom_list, flatten(atom_list)
                    else:
                        return atom_list, flatten(atom_list), annotate_info
                else:
                    if not annotate:
                        return atom_list
                    else:
                        return atom_list, annotate_info
        
        # Check if the number of atoms per step is correct
        if len(inp.split('_')) == 2 and inp.split('_')[0].lower() == "clear" and inp.split('_')[1].isdigit() and int(inp.split('_')[1]) < i:
            atom_list.pop(int(inp.split('_')[1])-1)
            print("Line %d is removed." % (int(inp.split('_')[1])))
            continue
        if len(inp.split()) + annotate != block_number:
            print(f"Warning: {block_number} integers seperated by {block_number-1} space(s).")
            continue
        
        if Atom_split_flag:
            tmp = []
        else:
            # Create an Atom object and set its bridge to the software object
            atoms = deepcopy(atom)
            atoms.bridge = software

        # Loop over the individual atoms in the input
        for a in inp.split()[:block_number]:
            # Check if the input is an integer
            if compile(r"^\-?\d+$").match(a) is None:
                print(f"Warning: Input {block_number} integers.")
                break
            if Atom_split_flag:
                # Create an Atom object and set its bridge to the software object
                atoms = deepcopy(atom)
                atoms.bridge = software
                err_return = atoms.put(a)
                # Check whether there is any error
                if atoms.err_list[err_return]:
                    print(atoms.err_list[err_return])
                else:
                    tmp.append(atoms)
        else:
            if Atom_split_flag:
                i += 1
                atom_list.append(tmp)
            else:
                # Get the return to check whether there is any error
                err_return = atoms.put(inp.split())

                # Check whether there is any error
                if atoms.err_list[err_return]:
                    print(atoms.err_list[err_return])
                else:
                    i += 1
                    atom_list.append(atoms)
        
        if annotate:
            annotate_info.append(inp.split()[block_number:block_number+annotate])

def blocks(software: object, block_number: int=2, items: int=-1, atom=Atom(), information: bool=True, flatten_flag: bool=False, annotate: int=0):
    """
    Prompt user to input lists of blocks.
    
    'blocks' function allows over one atom every block, and 'couple' function only allows an atom every block.
    
    if 'flatten_flag' is True, the function returns an additional 'Atom' object. The object contains without repetation.
    
    Args:
    -----
        software (mandatory): object, no default
            The software object with the necessary methods and attributes.
        
        block_number (required): int, default is 2
            The number of blocks to expect as input per item.
        
        items (required): int, default is -1
            The number of steps to expect as input. If -1, the input will continue until "end" is input.
        
        information (optional): bool, default is True
            Whether to display information about the start atom number.
        
        flatten_flag (optional): bool, default is False
            Whether to flatten the returned list.
    
    Returns:
    --------
        list: Return a two-dimensional list, including lists of items.
            Every item includes a list of blocks, whose number is equal to the 'block_number' parameter.
            Every block of inputted atoms uses an 'Atom' object, whose elements are sorted.
    
    Examples:
    ---------
    """
    
    # Initialize counter
    i = 1
    
    # Read the elements from the software object
    software.read_elements()
    atom_list = []
    
    # Read atom numbers from a file
    atom_list_from_file = file(software.args.atomfile, items)
    
    # Display the start atom number if information is True
    print(f"\nStart atom number: {software.args.atom}")
    if atom_list_from_file is None and information:
        print(f"Input {block_number} blocks separated by {block_number-1} space(s).\nInput numbers or elements separated by cammas in a region.\n")
    
    if annotate:
        annotate_info = []
    # Prompt the user to input blocks until the specified number of steps is reached or "end" is input
    while(1):
        if items > 0 and i > items:
            break
        
        if atom_list_from_file is None or atom_list_from_file == []:
            inp = input("%d: " % (i))
        else:
            inp = atom_list_from_file.pop(0)
            print(f"{i}: {inp}")
        
        if inp.lower() == "end":
            if items > -1:
                print("Warning: 'end' is invalid when the number of the step lines is specified.")
                continue
            else:
                if flatten_flag:
                    if annotate:
                        return atom_list, flatten(atom_list), annotate_info
                    else:
                        return atom_list, flatten(atom_list)
                else:
                    if annotate:
                        return atom_list, annotate_info
                    else:
                        return atom_list
        
        if len(inp.split('_')) == 2 and inp.split('_')[0].lower() == "clear" and inp.split('_')[1].isdigit() and int(inp.split('_')[1]) < i:
            blocks.pop(int(inp.split('_')[1])-1)
            print("Block %d is removed." % (int(inp.split('_')[1])))
            continue
        elif len(inp.split()) < block_number or len(inp.split()) > block_number + annotate:
            print(f"Warning: {block_number} blocks seperated by space(s).")
            continue
        
        blocks, recorded = [], []
        # Loop over the individual blocks in the input
        for block in inp.split()[:block_number]:

            # Create an Atom object and bridge from the software object
            atoms = deepcopy(atom)
            atoms.bridge = software
            
            # Loop over the individual atoms in each block
            flag = False
            for a in block.split(','):
                # Get the return to check whether there is any error
                if "rem" in a.lower():
                    for _ in range(len(atoms.elements)):
                        if _ not in recorded:
                            err_return = atoms.put(_)
                
                else:
                    err_return = atoms.put(a)
                recorded.extend(atoms.get())

                # Check whether there is any error
                if atoms.err_list[err_return]:
                    print(atoms.err_list[err_return])
                    flag = True
                    break
            else:
                atoms.sort()
                blocks.append(atoms)
                
            if flag:
                break
        else:
            if not atoms.err_list[err_return]:
                atom_list.append(tuple(blocks))
                i += 1
                if annotate:
                    annotate_info.append(inp.split()[block_number:block_number+annotate])

def line(software: object, items: int=-1, atom=Atom(), information: bool=True, annotate: int=0) -> Atom:
    """
    Prompt user to input a list of atoms. Return an 'Atom' object with inputted atom(s).
    
    Args:
    -----
        software (mandatory): object, no default
            The software object with the necessary methods and attributes.
        items (required): int, default is -1
            The number of steps to expect as input. If -1, the input will continue until "end" is input.
        information (optional): bool, default is True
            Whether to display information about the start atom number.
        
    Returns:
    --------
        Atom: An Atom object representing the input atoms.
    """
    
    # Initialize counter
    i = 1
    
    # Create an Atom object and bridge from the software object
    atoms = atom
    software.read_elements()
    atoms.bridge = software
    
    # Read atom numbers from a file
    atom_list_from_file = file(software.args.atomfile, items)
    
    # Display the start atom number if information is True
    print(f"\nStart atom number: {software.args.atom}")
    if atom_list_from_file is None and information:
        print("Input the numbers of atoms.\n")
    
    # Prompt the user to input atom numbers until the specified number of steps is reached or "end" is input
    annotate_info = []
    while(1):
        if items > 0 and i > items:
            break
        
        if atom_list_from_file is None or atom_list_from_file == []:
            inp = input("%d: " % (i))
        else:
            inp = atom_list_from_file.pop(0)
            print(f"{i}: {inp}")
        
        if inp.lower() == "end":
            if items > -1:
                print("Warning: 'end' is invalid when the number of the step lines is specified.")
                continue
            else:
                atoms.sort()
                return atoms
        
        if annotate:
            inp, tmp = inp.split()[0], inp.split()[1:1+annotate]
        # Get the return to check whether there is any error
        err_return = atoms.put(inp)
        
        # Check whether there is any error
        if atoms.err_list[err_return]:
            print(atoms.err_list[err_return])
        elif err_return == 2:
            if annotate:
                annotate_info.append(tmp)
                return atoms, annotate_info
            else:
                return atoms
        else:
            i += 1
            if annotate:
                annotate_info.append(tmp)

def lines(software, information: bool=True, items: int=-1, atom=Atom(), flatten_flag: bool=False, annotate: int=0):
    """
    Prompt the user to input a list of atom numbers. Return a list of 'Atom' objects. Every object includes an item.

    Parameters
    ----------
    software (mandatory): object, no default
        An object that has the attribute "preprocess" which is a dictionary.
    items (required): int, default is -1
        The number of steps to be input.
    information (optional): bool, default is True
        If True, print the information of start atom number.
    flatten_flag (optional): bool, default is False
        If True, return a flatten list of all atoms.

    Returns
    -------
    list or tuple
        A list of Atom objects representing the input atoms.
        If flatten_flag is True, return a tuple of two lists, the first list contains the atoms, 
        and the second list contains the flatten list of all atoms.
    """
    
    # Initialize counter
    i = 1  
    all_flag = False
    
    atom_list = []

    # Read atom numbers from a file
    atom_list_from_file = file(software.args.atomfile, items)
    
    # Print the start atom number if information is True
    print(f"\nStart atom number: {software.args.atom}")
    if atom_list_from_file is None and information:
        print("Input the numbers of atoms.\n")

    # Prompt the user to input atom numbers until the specified number of steps is reached or "end" is input
    annotate_info = []
    while(1):
        if items > 0 and i > items:
            # Add atoms without repetation to the flattened list if "flatten" is True
            if flatten_flag:
                if all_flag:
                    return (atom_list, "all", annotate_info) if annotate else (atom_list, "all")
                else:
                    return (atom_list, flatten(atom_list), annotate_info) if annotate else (atom_list, flatten(atom_list))
            else:
                return atom_list, annotate_info if annotate else atom_list
        
        if atom_list_from_file is None or atom_list_from_file == []:
            inp = input("%d: " % (i))
        else:
            inp = atom_list_from_file.pop(0)
            print(f"{i}: {inp}")

        if inp.lower() == "end":
            if items > -1:
                print("Warning: 'end' is invalid when the number of the step lines is specified.")
                continue
            else:
                if flatten_flag:
                    if all_flag:
                        return (atom_list, "all", annotate_info) if annotate else (atom_list, "all")
                    else:
                        return (atom_list, flatten(atom_list), annotate_info) if annotate else (atom_list, flatten(atom_list))
                else:
                    return (atom_list, annotate_info) if annotate else atom_list

        elif inp == "":
            print("Warning: Input error.")
            continue
        
        # Create an Atom object and bridge from the software object
        atoms = deepcopy(atom)
        atoms.bridge = software

        inp, tmp = inp.split()[0], inp.split()[1:1+annotate]
        if inp.lower() == "all":
            all_flag = True
            err_return = atoms.put(inp)
        elif inp.lower() == "same":
            err_return = atoms.put(atom_list[-1])
        elif len(inp.split('_')) == 2 and inp.split('_')[0].lower() == "clear" and inp.split('_')[1].isdigit() and int(inp.split('_')[1]) < i:
            atom_list.pop(int(inp.split('_')[1])-1)
            print("Line %d is removed." % (int(inp.split('_')[1])))
            continue
        else:
            sp = inp.split(',')
            for a in sp:
                err_return = atoms.put(a)
                if atoms.err_list[err_return]:
                    if atom_list_from_file is not None:
                        print(f"Line {i} in {software.args.atomfile}: {atoms.err_list[err_return]}")
                    else:
                        print(atoms.err_list[err_return])
                    break

        if not atoms.err_list[err_return]:
            atoms.sort()
            atom_list.append(atoms)
            annotate_info.append(tmp)
            i += 1

def file(filename, items=-1):
    if filename is not None:
        if path.isfile(filename):
            lines = []
            with open(filename) as read_file:
                idx = 1
                for line in read_file:
                    if line == '\n':
                        continue
                    idx += 1
                    lines.append(line.replace('\n', ''))
                    if items != -1 and idx > items:
                        break
            if items == -1 and lines[-1].lower() != "end":
                lines.append("end")
            return lines
        else:
            print(f"Warning: {filename} is not a file.")
            return None
    else:
        return None

def wrap(software, type="d"):
    if type == "d":
        from kit.function import Distance
        wrap = Distance(software)
    elif type == "a":
        from kit.function import Angle
        wrap = Angle(software)
    elif type == "dh":
        from kit.function import Dihedral_Angle
        wrap = Dihedral_Angle(software)
    else:
        raise NameError("'type' parameter should be 'd', 'a', or 'dh'.")
    wrap.bridge = software
    wrap.center_position = wrap.fractional_position[:, 0, :]
    wrap.measure_position = wrap.fractional_position[:, 1, :]
    
    while(1):
        yn = input("\nDoes measured atom wrap at step 0? (y/n) [y]: ")
        if yn == "" or yn.lower() == 'y':
            while(1):
                yn = input("Automatic, reversely automatic, sprinkler, or manual correction (a/ra/s/m) [a]: ")
                if yn == "" or yn.lower() == 'a':
                    wrap.auto()
                    break
                elif yn.lower() == 'ra':
                    wrap.rev_auto()
                    break
                elif yn.lower() == 's':
                    wrap.sprinkler()
                    break
                elif yn.lower() == 'm':
                    if isinstance(wrap, Distance):
                        wrap = manual_direction(1, wrap, 'd')
                    elif isinstance(wrap, Angle):
                        wrap = angle_manual(wrap)
                    elif isinstance(wrap, Dihedral_Angle):
                        wrap = dihedral_angle_manual(wrap)
                    wrap.manual()
                    break
                else:
                    print("Warning: Input 'a', 'ra', or 's'.\n")
        elif yn.lower() == 'n':
            wrap.center_position = wrap.wrap(wrap.center_position)
            wrap.measure_position = wrap.wrap(wrap.measure_position)
        else:
            print("Warning: Input 'y' or 'n'.\n")
            continue
        if type == "d":
            return wrap.distance()
        elif type == "a":
            return wrap.angle()
        elif type == "dh":
            return wrap.dihedral_angle()

def angle_manual(wrap):
    while(1):
        side_wrap = input("Modified side atom(s) (Over one atom allowed) (1/2): ")
        if (side_wrap != "1") ^ (side_wrap != "2") ^ (side_wrap != "1 2"):
            print("Warning: Input error.")
            continue
        if "1" in side_wrap.split():
            return manual_direction(1, wrap, 'a')
        if "2" in side_wrap.split():
            return manual_direction(2, wrap, 'a')

def dihedral_angle_manual(wrap):
    while(1):
        side_wrap = input("Modified side atom(s) (Over one atom allowed) (1/2/3): ")
        for num in side_wrap.split():
            if (num != "1") ^ (num != "2") ^ (num != "3"):
                print("Warning: Input error.")
                break
        if "1" in side_wrap.split():
            return manual_direction(1, wrap, 'a')
        if "2" in side_wrap.split():
            return manual_direction(2, wrap, 'a')
        if "3" in side_wrap.split():
            return manual_direction(3, wrap, 'a')

def manual_direction(No, software, kind):
    if kind == 'd':
        directions = input("Modified direction(s) (Over one direction allowed): ")
    else:
        directions = input(f"Modified direction(s) for m{No} (Over one direction allowed): ")
    for direction in directions.split():
        if 'a' in direction.lower():
            while(1):
                PNDirection = input("+/- A direction (+/-): ")
                if PNDirection == '+':
                    software.directions[No][0] += 1
                    break
                elif PNDirection == '-':
                    software.directions[No][0] -= 1
                    break
                else:
                    print("Warning: Input error.")
        if 'b' in direction.lower():
            while(1):
                PNDirection = input("+/- B direction (+/-): ")
                if PNDirection == '+':
                    software.directions[No][1] += 1
                    break
                elif PNDirection == '-':
                    software.directions[No][1] -= 1
                    break
                else:
                    print("Warning: Input error.")
        if 'c' in direction.lower():
            while(1):
                PNDirection = input("+/- C direction (+/-): ")
                if PNDirection == '+':
                    software.directions[No][2] += 1
                    break
                elif PNDirection == '-':
                    software.directions[No][2] -= 1
                    break
                else:
                    print("Warning: Input error.")
    return software

def threshold():
    while(1):
        yn = input("\nRevise the bond thresholds [n]: ")
        if yn == "":
            return None
        elif yn.lower() == 'y':
            break
        elif yn.lower() == 'n':
            return None
    if yn == 'y':
        i = 1
        bond_type = {}
        while(1):
            tmp = input(f"{i}: ")
            if tmp.lower() == "end":
                break
            sp = tmp.split()
            if len(sp) != 2:
                print("Warning: Input a bond type and the threshold separated by a space.")
                continue
            elif len(sp[0].split('-')) != 2:
                print("Warning: The format of bond type is two elements separated by '-', like 'C-C', 'C-O'.")
                continue
            elif compile(r"^\d+\.?\d*$").match(sp[1]) is None:
                print("Warning: The format of the threshold is a positive number.")
                continue
            i += 1
            bond_type[sp[0]] = float(sp[1])
            tmp = sp[0].split('-')
            bond_type[f"{tmp[1]}-{tmp[0]}"] = float(sp[1])
    return bond_type

def smooth(lines):
    from kit.function import Smooth
    from numpy import array
    while(1):
        yn = input("\nWould you like to smooth the curve (y/n): ")
        if yn.lower() == 'y':
            smooth = Smooth()
            while(1):
                tmp = input("Set the smooth factor: ")
                smooth.factor = tmp
                if smooth.factor != -1:
                    return array(smooth.smooth(lines, tmp)), smooth.factor
        elif yn.lower() == 'n':
            return lines, 0
        else:
            print("Warning: Input 'y' or 'n'.\n")
