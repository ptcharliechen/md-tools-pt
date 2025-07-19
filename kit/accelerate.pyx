import numpy as np
cimport numpy as np
from libc.math cimport sqrt

cpdef list EMA(np.ndarray unsmoothed_curves, float factor):
    cdef list curves = []
    cdef list smoothed_points = []
    cdef float former, point
    if unsmoothed_curves.ndim == 2:
        for unsmoothed_curve in unsmoothed_curves:
            smoothed_points = []
            for point in unsmoothed_curve:
                if smoothed_points:
                    former = smoothed_points[-1]
                    smoothed_points.append(former*factor + point*(1-factor))
                else:
                    smoothed_points.append(point)
            curves.append(smoothed_points)
        return curves
    elif unsmoothed_curves.ndim == 1:
        for point in unsmoothed_curves:
            if smoothed_points:
                former = smoothed_points[-1]
                smoothed_points.append(former*factor + point*(1-factor))
            else:
                smoothed_points.append(point)
        return smoothed_points

cpdef np.ndarray c2f_acc(np.ndarray cart_pos, np.ndarray lattice):
    cdef list frac_pos = []
    if lattice.ndim == 3:
        for pos, l in zip(cart_pos, lattice):
            frac_pos.append(c2f(pos, l))
    elif lattice.ndim == 2:
        for pos in cart_pos:
            frac_pos.append(c2f(pos, lattice))
    return np.array(frac_pos)

cpdef np.ndarray f2c_acc(np.ndarray frac_pos, np.ndarray lattice):
    cdef list cart_pos = []
    if lattice.ndim == 3:
        for pos, l in zip(frac_pos, lattice):
            cart_pos.append(f2c(pos, l))
    elif lattice.ndim == 2:
        for pos in frac_pos:
            cart_pos.append(f2c(pos, lattice))
    return np.array(cart_pos)

cpdef np.ndarray c2f(np.ndarray atom_pos, np.ndarray lattice):
    return atom_pos.dot(np.linalg.inv(lattice))

cpdef np.ndarray f2c(np.ndarray atom_pos, np.ndarray lattice):
    return atom_pos.dot(lattice)

cpdef np.ndarray shift_to_origin(np.ndarray atom_pos):
    cdef tuple indices

    while(np.any(atom_pos > 1) or np.any(atom_pos < 0)):
        indices = np.where(atom_pos > 1)
        if indices[0].size > 0:
            atom_pos[indices] -= 1

        indices = np.where(atom_pos < 0)
        if indices[0].size > 0:
            atom_pos[indices] += 1
    return atom_pos

cpdef np.ndarray image_shift(np.ndarray atom_pos, np.ndarray lattice, int benchmark=0, int cartesian=1):
    cdef int i, j

    if atom_pos.ndim == 2:
        for i in range(benchmark+1, len(atom_pos)):
            while(any((atom_pos[i] - atom_pos[i-1]) < -0.5)):
                idx = ((atom_pos[i] - atom_pos[i-1]) < -0.5)
                atom_pos[i][idx] += 1
            while(any((atom_pos[i] - atom_pos[i-1]) > 0.5)):
                idx = ((atom_pos[i] - atom_pos[i-1]) > 0.5)
                atom_pos[i][idx] -= 1
        if benchmark != 0:
            for i in range(benchmark-1, -1, -1):
                while(any((atom_pos[i] - atom_pos[i+1]) < -0.5)):
                    idx = ((atom_pos[i] - atom_pos[i+1]) < -0.5)
                    atom_pos[i][idx] += 1
                while(any((atom_pos[i] - atom_pos[i+1]) > 0.5)):
                    idx = ((atom_pos[i] - atom_pos[i+1]) > 0.5)
                    atom_pos[i][idx] -= 1
    elif atom_pos.ndim == 3:
        for i in range(1, len(atom_pos)):
            for j in range(len(atom_pos[0])):
                while(any((atom_pos[i][j] - atom_pos[i-1][j]) < -0.5)):
                    idx = ((atom_pos[i][j] - atom_pos[i-1][j]) < -0.5)
                    atom_pos[i][j][idx] += 1
                while(any((atom_pos[i][j] - atom_pos[i-1][j]) > 0.5)):
                    idx = ((atom_pos[i][j] - atom_pos[i-1][j]) > 0.5)
                    atom_pos[i][j][idx] -= 1

    if lattice.ndim == 2:
        if cartesian:
            return f2c(atom_pos, lattice)
        else:
            return atom_pos
    elif lattice.ndim == 3:
        if cartesian:
            return np.array(f2c_acc(atom_pos, lattice))
        else:
            return atom_pos

cpdef np.ndarray rev_image_shift(np.ndarray atom_pos, np.ndarray lattice):
    cdef int i, j

    atom_pos = shift_to_origin(atom_pos)

    if atom_pos.ndim == 2:
        for i in range(len(atom_pos)-2, -1, -1):
            while(any((atom_pos[i+1] - atom_pos[i]) < -0.5)):
                idx = ((atom_pos[i+1] - atom_pos[i]) < -0.5)
                atom_pos[i][idx] += 1
            while(any((atom_pos[i+1] - atom_pos[i]) > 0.5)):
                idx = ((atom_pos[i+1] - atom_pos[i]) > 0.5)
                atom_pos[i][idx] -= 1
    elif atom_pos.ndim == 3:
        for i in range(len(atom_pos[0])-2, -1, -1):
            for j in range(len(atom_pos)):
                while(any((atom_pos[i+1][j] - atom_pos[i][j]) < -0.5)):
                    idx = ((atom_pos[i+1][j] - atom_pos[i][j]) < -0.5)
                    atom_pos[i][j][idx] += 1
                while(any((atom_pos[i+1][j] - atom_pos[i][j]) > 0.5)):
                    idx = ((atom_pos[i+1][j] - atom_pos[i][j]) > 0.5)
                    atom_pos[i][j][idx] -= 1
    if lattice.ndim == 2:
        return f2c(atom_pos, lattice)
    elif lattice.ndim == 3:
        return np.array(f2c_acc(atom_pos, lattice))

cpdef np.ndarray angle(np.ndarray v1, np.ndarray v2, int degreeFlag=1):
    cdef int i
    cdef np.ndarray angles
    angles = np.zeros(len(v1))
    if degreeFlag:
        for i in range(len(v1)):
            angles[i] = np.degrees(np.arccos(v1[i].dot(v2[i])/(np.linalg.norm(v1[i])*np.linalg.norm(v2[i]))))
    else:
        for i in range(len(v1)):
            angles[i] = np.arccos(v1[i].dot(v2[i])/(np.linalg.norm(v1[i])*np.linalg.norm(v2[i])))
    return angles

cpdef np.ndarray dihedral_angle(np.ndarray v1, np.ndarray v2, np.ndarray v3, int degreeFlag=1):
    cdef int i
    cdef float tmp
    cdef np.ndarray dh_angles, n1, n2
    n1 = np.cross(v1, v2)
    n2 = np.cross(v1, v3)
    dh_angles = np.zeros(len(n1))
    if degreeFlag:
        for i in range(len(n1)):
            tmp = np.degrees(np.arccos(n1[i].dot(n2[i])/(np.linalg.norm(n1[i])*np.linalg.norm(n2[i]))))
            dh_angles[i] = (tmp if tmp < 90 else 180 - tmp)
    else:
        for i in range(len(n1)):
            tmp = np.arccos(n1[i].dot(n2[i])/(np.linalg.norm(n1[i])*np.linalg.norm(n2[i])))
            dh_angles[i] = (tmp if tmp < np.pi/2 else np.pi - tmp)
    return dh_angles

cpdef np.ndarray distance_matrix_wrap(np.ndarray cen_pos, np.ndarray mea_pos, np.ndarray lattice, np.ndarray wrap):
    cdef np.ndarray[np.float64_t, ndim=4] cen, mea
    cdef np.ndarray[np.int_t, ndim=4] wrap_arr
    cdef np.ndarray[np.float64_t, ndim=4] vectors

    cen_pos = cen_pos[:, None, None, :]  # (N_cen, 1, 1, 3)
    mea_pos = mea_pos[None, :, None, :]  # (1, N_mea, 1, 3)
    wrap_arr = wrap[None, None, :, :]        # (1, 1, N_wrap, 3)

    cen, mea = grid(cen_pos, mea_pos, lattice, wrap_arr)

    vectors = cen - mea
    return np.linalg.norm(vectors, axis=3)


cpdef tuple grid(np.ndarray cen_pos, np.ndarray mea_pos, np.ndarray lattice, np.ndarray wrap):
    cdef np.ndarray[np.float64_t, ndim=4] cen_pos_new, mea_pos_new

    cen_pos_new = np.matmul(cen_pos, lattice)           # (N_cen, 1, 1, 3) @ (3, 3) → (N_cen, 1, 1, 3)
    mea_pos_new = np.matmul(mea_pos + wrap, lattice)    # (1, N_mea, N_wrap, 3) @ (3, 3) → (1, N_mea, N_wrap, 3)

    return (cen_pos_new, mea_pos_new)

cpdef np.ndarray distance_matrix(np.ndarray[np.float64_t, ndim=2] cen_pos,
                                 np.ndarray[np.float64_t, ndim=2] mea_pos,
                                 np.ndarray[np.float64_t, ndim=2] lattice):
    cdef int i, j, m, n
    cdef int cen_shape = cen_pos.shape[0]
    cdef int mea_shape = mea_pos.shape[0]

    # memoryview
    cdef double[:, :] cen_mv = cen_pos
    cdef double[:, :] mea_mv = mea_pos
    cdef double[:, :] lat_mv = lattice

    # output
    cdef np.ndarray[np.float64_t, ndim=2] distances = np.empty((cen_shape, mea_shape), dtype=np.float64)
    cdef double[:, :] dist_mv = distances

    # intermediate vector
    cdef double cen_latt[3], mea_latt[3], diff[3]
    cdef double dx, dy, dz

    for i in range(cen_shape):
        for m in range(3):
            cen_latt[m] = 0.0
            for n in range(3):
                cen_latt[m] += cen_mv[i, n] * lat_mv[n, m]

        for j in range(mea_shape):
            for m in range(3):
                mea_latt[m] = 0.0
                for n in range(3):
                    mea_latt[m] += mea_mv[j, n] * lat_mv[n, m]

                diff[m] = cen_latt[m] - mea_latt[m]

            dx, dy, dz = diff[0], diff[1], diff[2]
            dist_mv[i, j] = sqrt(dx*dx + dy*dy + dz*dz)

    return distances
