import numpy as np
cimport numpy as cnp
cimport cython

cnp.import_array()

FTYPE = np.float64
CTYPE = np.complex128

ctypedef cnp.float64_t FTYPE_t
ctypedef cnp.complex128_t CTYPE_t


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def _compute_r_term(
    int m,
    cnp.ndarray[FTYPE_t, ndim=1] k,
    int ns,
    cnp.ndarray[CTYPE_t, ndim=2] Q,
    cnp.ndarray[CTYPE_t, ndim=2] brt,
    cnp.ndarray[FTYPE_t, ndim=1] lam,
    cnp.ndarray[FTYPE_t, ndim=1] ffm,
    int nr,
    cnp.ndarray[FTYPE_t, ndim=1] ffp,
    cnp.ndarray[CTYPE_t, ndim=3] psi,
    cnp.ndarray[CTYPE_t, ndim=2] psir,
    float rss,
    cnp.ndarray brt_outer,
):
    for l in range(ns):
        # Ignore the l=0 and m=0 term; for a globally divergence free field
        # this term is zero anyway, but numerically it may be small which
        # causes numerical issues when solving for c, d
        if l == 0 and m == 0:
            continue
        # - sum (c_{lm} + d_{lm}) * lam_{l}
        # lam[l] is small so this blows up
        cdlm = np.dot(Q[:, l], brt[:, m]) / lam[l]

        if brt_outer is None:
            # - ratio c_{lm}/d_{lm} [numerically safer this way up]
            ratio = (ffm[l] ** (nr - 1) - ffm[l] ** nr) / (ffp[l] ** nr - ffp[l] ** (nr - 1))
            dlm = cdlm / (1.0 + ratio)
            clm = ratio * dlm
        else:
            cdlm1 = np.dot(Q[:, l], brt_outer[:, m]) / lam[l] * rss**2
            clm = (cdlm1 - ffm[l] ** nr * cdlm) / (ffp[l] ** nr - ffm[l] ** nr)
            dlm = (cdlm1 - ffp[l] ** nr * cdlm) / (ffm[l] ** nr - ffp[l] ** nr)

        psir[:, l] = clm * ffp[l] ** k + dlm * ffm[l] ** k

    # - compute entry for this m in psit = Sum_l c_{lm}Q_{lm}**j
    psi[:, :, m] = np.dot(psir, Q.T)
    return psi, psir


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def _als_alp(
    int nr,
    int nphi,
    cnp.ndarray[FTYPE_t, ndim=1] Fs,
    cnp.ndarray[FTYPE_t, ndim=3] psi,
    cnp.ndarray[FTYPE_t, ndim=1] Fp,
    cnp.ndarray[FTYPE_t, ndim=3] als,
    cnp.ndarray[FTYPE_t, ndim=3] alp,
):
    for j in range(nr + 1):
        for i in range(nphi + 1):
            als[i, :, j] = Fs * (psi[j, :, ((i - 1) % nphi)] - psi[j, :, ((i) % nphi)])
        for i in range(nphi):
            alp[i, 1:-1, j] = Fp[1:-1] * (psi[j, 1:, i] - psi[j, :-1, i])
    return als, alp


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def _A_diag(
    cnp.ndarray[FTYPE_t, ndim=2] A,
    int ns,
    cnp.ndarray[FTYPE_t, ndim=1] Vg,
    cnp.ndarray[FTYPE_t, ndim=1] Uc,
    cnp.ndarray[FTYPE_t, ndim=1] mu,
    int m,
):
    for j in range(ns):
        A[j, j] = Vg[j] + Vg[j + 1] + Uc[j] * mu[m]
    return A
