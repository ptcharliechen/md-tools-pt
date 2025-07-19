from os import path, getcwd
from numpy import zeros, array, arange, ndarray, float64

def read_csv(output, index_flag=True, column_flag=True, num_flag=False):
    data = []
    if index_flag:
        index = []
    else:
        index = None
    if not column_flag:
        column = None
    with open(output) as read_file:
        for row, line in enumerate(read_file):
            if row == 0 and column_flag:
                column = line.replace('\n', '').split(',')
            else:
                sp = line.replace('\n', '').split(',')
                if index_flag:
                    if sp[0].isdigit():
                        index.append(int(sp[0]))
                    else:
                        index.append(sp[0])
                    data.append(sp[1:])
                else:
                    data.append(sp)
        if num_flag:
            numData = []
            for row_data in data:
                numData.append([float(num) for num in row_data])
            return numData, index, column
        else:
            return data, index, column

def write_csv(output, data, indices=None, columns=None):
    from collections.abc import Iterable
    if path.splitext(output)[1] != ".csv":
        output += ".csv"
    with open(path.join(getcwd(), output), "w") as write_file:
        if isinstance(columns, Iterable) and not isinstance(columns, str):
            for column in columns:
                write_file.write(",{}".format(column))
            write_file.write("\n")
        elif isinstance(columns, str):
            write_file.write(",{}\n".format(columns))
        if isinstance(columns, Iterable) and isinstance(data, ndarray):
            if data.ndim == 2 and data.shape[1] > len(columns) * 30:
                data = data.T
        if indices is not None:
            for idx, datum in enumerate(data):
                write_file.write("{}".format(indices[idx]))
                if isinstance(datum, ndarray):
                    for c in datum:
                        write_file.write(",{:.4f}".format(c))
                elif isinstance(datum, float64) or isinstance(datum, float):
                    write_file.write(",{:.4f}".format(datum))
                elif isinstance(datum, Iterable):
                    if isinstance(datum[0], Iterable) and not isinstance(datum[0], str):
                        for c in datum:
                            write_file.write(",{}".format(c))
                    else:
                        write_file.write(",{}".format(datum))
                write_file.write('\n')
        else:
            for datum in data:
                if isinstance(datum, ndarray):
                    for c in datum:
                        write_file.write(",{:.4f}".format(c))
                elif isinstance(datum, float64) or isinstance(datum, float64):
                    write_file.write(",{:.4f}".format(datum))
                elif isinstance(datum, Iterable):
                    if isinstance(datum[0], Iterable) and not isinstance(datum[0], str):
                        for c in datum:
                            write_file.write(",{}".format(c))
                    else:
                        write_file.write(",{}".format(datum))
                write_file.write('\n')

def cubic_spline(x, y, decimals=4, below_zero=True):
    h = zeros(len(x)-1)
    for i in range(len(h)):
        h[i] = x[i+1] - x[i]
    A = zeros((len(x), len(x)))
    r = zeros(len(x))
    A[-1, -1] = 1
    A[0, 0] = 1
    
    for i in range(1, len(A)-1):
        A[i, i] = 2*(h[i-1] + h[i])
        A[i, i-1] = h[i-1]
        A[i, i+1] = h[i]
        r[i] = 3/h[i] * (y[i+1]-y[i]) - 3/h[i-1] * (y[i]-y[i-1])

    coef = zeros((len(x), 4))
    coef[:, 0] = array(y).T
    coef[:, 2] = Thomas(A, r)
    for i in range(len(x)-1):
        coef[i, 1] = (coef[i+1, 0]-coef[i, 0])/h[i] - 2*coef[i, 2]*h[i]/3 - coef[i+1, 2]*h[i]/3
        coef[i, 3] = (coef[i+1, 2]-coef[i, 2])/(3*h[i])

    x_val = arange(max(x))
    y_val = zeros(len(x_val))
    for i in range(len(x)-1):
        for j in range(int(x[i])-int(x[0]), int(x[i+1])-int(x[0])):
            y_val[j] = coef[i, 0] + coef[i, 1]*(x_val[j]-x[i]) + coef[i, 2]*(x_val[j]-x[i])**2 + coef[i, 3]*(x_val[j]-x[i])**3
            if not below_zero and y_val[j] < 0:
                y_val[j] = 0

    if decimals > 0:
        from numpy import round
        y_val = round(y_val, decimals)
        y_val[-1] = round(y[-1], decimals)
    else:
        y_val[-1] = y[-1]

    return y_val

def Thomas(A, r):
    x = zeros(A.shape[0])
    if abs(A[0, 1]) < 10**(-8):
        start = 1
        r[1] -= r[0]*A[1, 0]
        x[0] = r[0]/A[0, 0]
    else:
        start = 0
    end = A.shape[0]
    
    for i in range(start+1, end):
        A[i, i-1] /= A[i-1, i-1]
        A[i, i] -= A[i, i-1]*A[i-1, i]
        r[i] -= A[i, i-1]*r[i-1]
    
    x[A.shape[0]-1] = r[A.shape[0]-1]/A[A.shape[0]-1, A.shape[0]-1]
    for i in range(end-2, start-1, -1):
        x[i] = (r[i] - A[i, i+1]*x[i+1])/A[i, i]
    
    return x
