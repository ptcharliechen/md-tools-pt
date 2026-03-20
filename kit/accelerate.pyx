import numpy as np
cimport numpy as cnp
cimport cython
from libc.math cimport sqrt, INFINITY, NAN
from cython.parallel cimport prange

cpdef list EMA(cnp.ndarray unsmoothed_curves, float factor):
    cdef list curves = []
    cdef list smoothed_points = []
    cdef float former, point
    if unsmoothed_curves.ndim == 2:
        unsmoothed_curves = unsmoothed_curves.T
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

cpdef cnp.ndarray c2f_acc(cnp.ndarray cart_pos, cnp.ndarray lattice):
    cdef list frac_pos = []
    if lattice.ndim == 3:
        for pos, l in zip(cart_pos, lattice):
            frac_pos.append(c2f(pos, l))
    elif lattice.ndim == 2:
        for pos in cart_pos:
            frac_pos.append(c2f(pos, lattice))
    return np.array(frac_pos)

cpdef cnp.ndarray f2c_acc(cnp.ndarray frac_pos, cnp.ndarray lattice):
    cdef list cart_pos = []
    if lattice.ndim == 3:
        for pos, l in zip(frac_pos, lattice):
            cart_pos.append(f2c(pos, l))
    elif lattice.ndim == 2:
        for pos in frac_pos:
            cart_pos.append(f2c(pos, lattice))
    return np.array(cart_pos)

cpdef cnp.ndarray c2f(cnp.ndarray atom_pos, cnp.ndarray lattice):
    return atom_pos.dot(np.linalg.inv(lattice))

cpdef cnp.ndarray f2c(cnp.ndarray atom_pos, cnp.ndarray lattice):
    return atom_pos.dot(lattice)

cpdef cnp.ndarray shift_to_origin(cnp.ndarray atom_pos):
    cdef tuple indices

    while(np.any(atom_pos > 1) or np.any(atom_pos < 0)):
        indices = np.where(atom_pos > 1)
        if indices[0].size > 0:
            atom_pos[indices] -= 1

        indices = np.where(atom_pos < 0)
        if indices[0].size > 0:
            atom_pos[indices] += 1
    return atom_pos

cpdef cnp.ndarray image_shift(cnp.ndarray atom_pos, cnp.ndarray lattice, int benchmark=0, int cartesian=1):
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
    
    if cartesian:
        if lattice.ndim == 2:
            return f2c(atom_pos, lattice)
        elif lattice.ndim == 3:
            return np.array(f2c_acc(atom_pos, lattice))
    else:
        return atom_pos

cpdef cnp.ndarray rev_image_shift(cnp.ndarray atom_pos, cnp.ndarray lattice):
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

cpdef cnp.ndarray angle(cnp.ndarray v1, cnp.ndarray v2, bint degree=True):
    cdef int i
    cdef cnp.ndarray angles
    angles = np.zeros(len(v1))
    if degree:
        for i in range(len(v1)):
            angles[i] = np.degrees(np.arccos(v1[i].dot(v2[i])/(np.linalg.norm(v1[i])*np.linalg.norm(v2[i]))))
    else:
        for i in range(len(v1)):
            angles[i] = np.arccos(v1[i].dot(v2[i])/(np.linalg.norm(v1[i])*np.linalg.norm(v2[i])))
    return angles

cpdef cnp.ndarray dihedral_angle(cnp.ndarray v1, cnp.ndarray v2, cnp.ndarray v3, bint degree=True):
    cdef int i
    cdef float angle
    cdef cnp.ndarray dh_angles, n1, n2

    n1 = np.cross(v1, v2)
    n2 = np.cross(v1, v3)
    dh_angles = np.zeros(len(n1))
    if degree:
        for i in range(len(n1)):
            angle = np.degrees(np.arccos(n1[i].dot(n2[i])/(np.linalg.norm(n1[i])*np.linalg.norm(n2[i]))))
            dh_angles[i] = (angle if angle < 90 else 180 - angle)
    else:
        for i in range(len(n1)):
            angle = np.arccos(n1[i].dot(n2[i])/(np.linalg.norm(n1[i])*np.linalg.norm(n2[i])))
            dh_angles[i] = (angle if angle < np.pi/2 else np.pi - angle)
    return dh_angles

#cpdef cnp.ndarray distance_matrix_wrap(cnp.ndarray cen_pos, cnp.ndarray mea_pos, cnp.ndarray lattice, cnp.ndarray wrap):
#    cdef cnp.ndarray[cnp.float64_t, ndim=4] cen, mea
#    cdef cnp.ndarray[cnp.int_t, ndim=4] wrap_arr
#    cdef cnp.ndarray[cnp.float64_t, ndim=4] vectors
#
#    cen_pos = cen_pos[:, None, None, :]
#    mea_pos = mea_pos[None, :, None, :]
#    wrap_arr = wrap[None, None, :, :]
#
#    cen, mea = grid(cen_pos, mea_pos, lattice, wrap_arr)
#
#    vectors = cen - mea
#    return np.linalg.norm(vectors, axis=3)

cpdef cnp.ndarray distance_matrix_wrap(cnp.ndarray cen_pos, cnp.ndarray mea_pos, cnp.ndarray lattice, cnp.ndarray wrap):
    cdef cnp.ndarray[cnp.float64_t, ndim=4] cen, mea
    cdef cnp.ndarray[cnp.int_t, ndim=4] wrap_arr
    cdef cnp.ndarray[cnp.float64_t, ndim=4] vectors

    # 使用 broadcasting 擴展 cen_pos 和 mea_pos 的 shape
    cen_pos = cen_pos[:, None, None, :]  # (N_cen, 1, 1, 3)
    mea_pos = mea_pos[None, :, None, :]  # (1, N_mea, 1, 3)
    wrap_arr = wrap[None, None, :, :]        # (1, 1, N_wrap, 3)

    # 傳給 grid（cen_pos 和 mea_pos 已 broadcast 好）
    cen, mea = grid(cen_pos, mea_pos, lattice, wrap_arr)

    vectors = cen - mea
    return np.linalg.norm(vectors, axis=3)

cpdef tuple grid(cnp.ndarray cen_pos, cnp.ndarray mea_pos, cnp.ndarray lattice, cnp.ndarray wrap):
    cdef cnp.ndarray[cnp.float64_t, ndim=4] cen_pos_new, mea_pos_new

    # 直接使用 broadcasting 加總後乘 lattice
    cen_pos_new = np.matmul(cen_pos, lattice)           # (N_cen, 1, 1, 3) @ (3, 3) → (N_cen, 1, 1, 3)
    mea_pos_new = np.matmul(mea_pos + wrap, lattice)    # (1, N_mea, N_wrap, 3) @ (3, 3) → (1, N_mea, N_wrap, 3)

    return (cen_pos_new, mea_pos_new)

#cdef void grid(double[:, :, :, :] cen_arr,
#                       double[:, :, :, :] mea_arr,
#                       long[:, :, :, :] wrap_arr,
#                       double[:, :] lattice,
#                       double[:, :, :, :] cen_out,
#                       double[:, :, :, :] mea_out):
#    cdef int i, j, k, m, n
#    cdef double shift[3]
#
#    for i in range(cen_arr.shape[0]):
#        for j in range(mea_arr.shape[1]):
#            for k in range(wrap_arr.shape[2]):
#                for n in range(3):
#                    # wrap shift = lattice.T @ wrap
#                    shift[n] = 0.0
#                    for m in range(3):
#                        shift[n] += lattice[m, n] * wrap_arr[0, 0, k, m]
#
#                    # cen = cen_arr @ lattice
#                    cen_out[i, j, k, n] = 0.0
#                    mea_out[i, j, k, n] = shift[n]
#                    for m in range(3):
#                        cen_out[i, j, k, n] += cen_arr[i, 0, 0, m] * lattice[m, n]
#                        mea_out[i, j, k, n] += mea_arr[0, j, 0, m] * lattice[m, n]
#
#
#cpdef cnp.ndarray distance_matrix_wrap(cnp.ndarray[cnp.float64_t, ndim=2] cen_pos,
#                                      cnp.ndarray[cnp.float64_t, ndim=2] mea_pos,
#                                      cnp.ndarray[cnp.float64_t, ndim=2] lattice,
#                                      cnp.ndarray[cnp.int64_t, ndim=2] wrap):
#    cdef int cen_shape = cen_pos.shape[0]
#    cdef int mea_shape = mea_pos.shape[0]
#    cdef int wrap_len  = wrap.shape[0]
#    cdef int i, j, k, d
#
#    # memoryview for lattice
#    cdef double[:, :] lattice_mv = lattice
#
#    # broadcasted inputs (manually simulated)
#    cdef cnp.ndarray[cnp.float64_t, ndim=4] cen_arr = np.empty((cen_shape, 1, 1, 3), dtype=np.float64)
#    cdef cnp.ndarray[cnp.float64_t, ndim=4] mea_arr = np.empty((1, mea_shape, 1, 3), dtype=np.float64)
#    cdef cnp.ndarray[cnp.int64_t, ndim=4] wrap_arr = np.empty((1, 1, wrap_len, 3), dtype=np.int64)
#    cen_arr[:, 0, 0, :] = cen_pos
#    mea_arr[0, :, 0, :] = mea_pos
#    wrap_arr[0, 0, :, :] = wrap
#
#    # output arrays for transformed positions
#    cdef cnp.ndarray[cnp.float64_t, ndim=4] cen = np.empty((cen_shape, mea_shape, wrap_len, 3), dtype=np.float64)
#    cdef cnp.ndarray[cnp.float64_t, ndim=4] mea = np.empty((cen_shape, mea_shape, wrap_len, 3), dtype=np.float64)
#
#    grid(cen_arr, mea_arr, wrap_arr, lattice_mv, cen, mea)
#
#    # compute distances
#    cdef cnp.ndarray[cnp.float64_t, ndim=3] dists = np.empty((cen_shape, mea_shape, wrap_len), dtype=np.float64)
#    cdef double[:, :, :, :] cen_mv = cen
#    cdef double[:, :, :, :] mea_mv = mea
#    cdef double[:, :, :] dists_mv = dists
#
#    cdef double dx, dy, dz
#
#    for i in range(cen_shape):
#        for j in range(mea_shape):
#            for k in range(wrap_len):
#                dx = cen_mv[i, j, k, 0] - mea_mv[i, j, k, 0]
#                dy = cen_mv[i, j, k, 1] - mea_mv[i, j, k, 1]
#                dz = cen_mv[i, j, k, 2] - mea_mv[i, j, k, 2]
#                dists_mv[i, j, k] = sqrt(dx*dx + dy*dy + dz*dz)
#
#    return dists

cpdef cnp.ndarray distance_matrix(cnp.ndarray[cnp.float64_t, ndim=2] cen_pos,
                                 cnp.ndarray[cnp.float64_t, ndim=2] mea_pos,
                                 cnp.ndarray[cnp.float64_t, ndim=2] lattice):
    cdef int i, j, m, n
    cdef int cen_shape = cen_pos.shape[0]
    cdef int mea_shape = mea_pos.shape[0]

    # memoryview
    cdef double[:, :] cen_mv = cen_pos
    cdef double[:, :] mea_mv = mea_pos
    cdef double[:, :] lat_mv = lattice

    # output
    cdef cnp.ndarray[cnp.float64_t, ndim=2] distances = np.empty((cen_shape, mea_shape), dtype=np.float64)
    cdef double[:, :] dist_mv = distances

    # intermediate vector
    cdef double cen_latt[3]
    cdef double mea_latt[3]
    cdef double diff[3]
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
#cpdef cnp.ndarray distance_matrix(cnp.ndarray cen_pos, cnp.ndarray mea_pos, cnp.ndarray lattice):
#    cdef cnp.ndarray[np.float64_t, ndim=3] cen, mea
#    cdef cnp.ndarray[np.float64_t, ndim=3] vectors
#
#    # Broadcasting 展開維度：
#    cen_pos = cen_pos[:, None, :]  # (N_cen, 1, 3)
#    mea_pos = mea_pos[None, :, :]  # (1, N_mea, 3)
#
#    # 將兩者乘上 lattice，得到在 lattice 空間中的座標
#    cen = np.matmul(cen_pos, lattice)  # (N_cen, 1, 3) → (N_cen, 1, 3)
#    mea = np.matmul(mea_pos, lattice)  # (1, N_mea, 3) → (1, N_mea, 3)
#
#    # 計算差向量與距離
#    vectors = cen - mea  # shape (N_cen, N_mea, 3)
#    return np.linalg.norm(vectors, axis=2)

ctypedef fused real_t:
    cnp.float64_t
    cnp.float32_t

@cython.cfunc
@cython.inline
cdef bint _axis_allows(cnp.int_t t, double rc, double L, double cutoff) nogil:
    if t == 0:
        return True
    if t == -1:
        return rc < cutoff
    return (L - rc) < cutoff  # t == +1

@cython.cfunc
@cython.inline
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cdef bint _image_allowed(const cnp.int_t[:] t, const double[:] rc,
                         const double[:] L, const double cutoff) nogil:
    return (_axis_allows(t[0], rc[0], L[0], cutoff) and
            _axis_allows(t[1], rc[1], L[1], cutoff) and
            _axis_allows(t[2], rc[2], L[2], cutoff))


cpdef tuple distance_matrix_cutoff(
    cnp.ndarray[real_t, ndim=2] center_positions, 
    cnp.ndarray[real_t, ndim=2] measure_positions, 
    cnp.ndarray[real_t, ndim=2] lattice, 
    cnp.ndarray[cnp.int_t, ndim=2] images, 
    double cutoff,
    bint minimal_only,
    int nthreads=4
):
    cdef Py_ssize_t Nc=center_positions.shape[0], Nm=measure_positions.shape[0], M=images.shape[0]
    center_positions = np.ascontiguousarray(center_positions)
    measure_positions = np.ascontiguousarray(measure_positions)
    lattice = np.ascontiguousarray(lattice)
    images = np.ascontiguousarray(images, dtype=np.int64)

    cdef real_t[:, :] cen = center_positions
    cdef real_t[:, :] mea = measure_positions
    cdef real_t[:, :] lat = lattice
    cdef cnp.int_t[:, :] imgs = images

    # Lx, Ly, Lz for directional cutoff
    cdef double[:] L = np.empty(3, dtype=np.float64)
    L[0] = <double>lat[0,0]
    L[1] = <double>lat[1,1]
    L[2] = <double>lat[2,2]

    cdef Py_ssize_t i, j, m, best_m
    cdef double[:] rc = np.empty(3, dtype=np.float64)
    cdef double dx, dy, dz, pos_x, pos_y, pos_z, d2, best_d2
    cdef bint any_allowed
    cdef cnp.ndarray npD
    cdef double[:, :] D2
    cdef double[:, :, :] D3
    cdef cnp.ndarray npT
    cdef cnp.int_t[:, :] T2

    if minimal_only:
        npD = np.empty((Nc, Nm), dtype=np.float64)
        npT = np.zeros((Nc, Nm), dtype=np.int64)
        D2 = npD
        T2 = npT

        with nogil:
            for i in prange(Nc, schedule='static', num_threads=nthreads):
                rc[0] = lat[0, 0]*<double>cen[i,0] + lat[0, 1]*<double>cen[i,1] + lat[0, 2]*<double>cen[i,2]
                rc[1] = lat[1, 0]*<double>cen[i,0] + lat[1, 1]*<double>cen[i,1] + lat[1, 2]*<double>cen[i,2]
                rc[2] = lat[2, 0]*<double>cen[i,0] + lat[2, 1]*<double>cen[i,1] + lat[2, 2]*<double>cen[i,2]
                for j in range(Nm):
                    best_d2 = INFINITY; best_m = -1; any_allowed = False
                    for m in range(M):
                        if _image_allowed(imgs[m], rc, L, cutoff):
                            any_allowed = True
                            dx = <double>mea[j,0] + imgs[m,0] - <double>cen[i,0]
                            dy = <double>mea[j,1] + imgs[m,1] - <double>cen[i,1]
                            dz = <double>mea[j,2] + imgs[m,2] - <double>cen[i,2]
                            pos_x = lat[0, 0]*dx + lat[0, 1]*dy + lat[0, 2]*dz
                            pos_y = lat[1, 0]*dx + lat[1, 1]*dy + lat[1, 2]*dz
                            pos_z = lat[2, 0]*dx + lat[2, 1]*dy + lat[2, 2]*dz
                            d2 = pos_x*pos_x + pos_y*pos_y + pos_z*pos_z
                            if d2 < best_d2:
                                best_d2 = d2; best_m = m
                    T2[i, j] = best_m
                    if any_allowed and best_m >= 0:
                        D2[i, j] = sqrt(best_d2)
                        T2[i, j] = best_m
                    else:
                        D2[i, j] = NAN

        return (npD, npT)

    else:
        npD = np.empty((Nc, Nm, M), dtype=np.float64)
        npD.fill(np.nan)
        D3 = npD
        npT = np.array(images, copy=True, dtype=np.int64)

        with nogil:
            for i in prange(Nc, schedule='static', num_threads=nthreads):
                rc[0] = lat[0, 0]*<double>cen[i, 0] + lat[0, 1]*<double>cen[i, 1] + lat[0, 2]*<double>cen[i, 2]
                rc[1] = lat[1, 0]*<double>cen[i, 0] + lat[1, 1]*<double>cen[i, 1] + lat[1, 2]*<double>cen[i, 2]
                rc[2] = lat[2, 0]*<double>cen[i, 0] + lat[2, 1]*<double>cen[i, 1] + lat[2, 2]*<double>cen[i, 2]
                for j in range(Nm):
                    for m in range(M):
                        if _image_allowed(imgs[m], rc, L, cutoff):
                            dx = <double>mea[j,0] + imgs[m,0] - <double>cen[i,0]
                            dy = <double>mea[j,1] + imgs[m,1] - <double>cen[i,1]
                            dz = <double>mea[j,2] + imgs[m,2] - <double>cen[i,2]
                            pos_x = lat[0, 0]*dx + lat[0, 1]*dy + lat[0, 2]*dz
                            pos_y = lat[1, 0]*dx + lat[1, 1]*dy + lat[1, 2]*dz
                            pos_z = lat[2, 0]*dx + lat[2, 1]*dy + lat[2, 2]*dz
                            D3[i, j, m] = sqrt(pos_x*pos_x + pos_y*pos_y + pos_z*pos_z)
        return (npD, np.nan)
