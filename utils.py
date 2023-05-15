"""
Helper functions for the implementation of the Riemannian logarithm on flag manifolds.

All the details can be found in [Szwagier2023].

Lead author: Tom Szwagier.

References
----------
.. [Szwagier2023] T. Szwagier, X. Pennec.
    “Rethinking the Riemannian logarithm on flag manifolds as an orthogonal alignment problem.”
    To appear: GSI'23 6th International Conference on Geometric Science of Information.
"""

from geomstats.geometry.skew_symmetric_matrices import SkewSymmetricMatrices
import itertools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ortho_group, special_ortho_group
from scipy.linalg import expm, logm, block_diag, svd
from time import time


def projector_embedding(Q, signature_0):
    r""" Embed a flag in a product of Grassmannians, by forming orthogonal projection
    matrices relative to each subspace of the flag. Corresponds to the embedding map
    :math:`\Pi_{\operatorname{Gr}(\I)}` in [Szwagier2023].

    Parameters
    ----------
    Q: array-like, shape=[n, n]
        An orthogonal matrix representing a flag.
    signature_0: tuple of ints
        The signature of the flag.

    Returns
    -------
    Q_proj: array-like, shape=[n, n * r]
        The sequence of orthogonal projection matrices related to the flag.
        The sequence is horizontally concatenated into one (n, n * r) array.
    """
    n = signature_0[-1]
    Q_proj = np.zeros((signature_0[-1], signature_0[-1] * len(signature_0[1:])))
    for i in range(len(signature_0[1:])):
        Q_i = Q[:, signature_0[i]: signature_0[i + 1]]
        Q_proj[:, n * i:n * (i + 1)] = (Q_i @ Q_i.T)
    return Q_proj


def projection_kernel(signature_0):
    """ Create a binary mask with 1s on the diagonal blocks related to a given
    signature.

    Parameters
    ----------
    signature_0: tuple of ints
        The signature of the flag.

    Returns
    -------
    skew_diag_kernel: array-like, shape=[n, n]
        A binary mask with 1s on the diagonal blocks related to a given signature.
    """
    n = signature_0[-1]
    skew_diag_kernel = np.zeros((n, n))
    for (d, d_) in zip(signature_0[:-1], signature_0[1:]):
        skew_diag_kernel[d:d_, d:d_] = np.ones((d_ - d, d_ - d))
    return skew_diag_kernel


def proj_H(X, skew_diag_kernel):
    r""" Project a tangent vector of :math:`\mathcal{O}(n)` onto the horizontal space
    w.r.t the flag action.

    Parameters
    ----------
    X: array-like, shape=[n, n]
        A tangent vector in O(n).
    skew_diag_kernel: array-like, shape=[n, n]
        A binary mask with 1s on the diagonal blocks related to a given signature.

    Returns
    -------
    H: array-like, shape=[n, n]
        A horizontal vector.
    """
    return X * (1 - skew_diag_kernel)


def proj_V(X, skew_diag_kernel):
    r""" Project a tangent vector of :math:`\mathcal{O}(n)` onto the vertical space
    w.r.t the flag action.

    Parameters
    ----------
    X: array-like, shape=[n, n]
        A tangent vector in O(n).
    skew_diag_kernel: array-like, shape=[n, n]
        A binary mask with 1s on the diagonal blocks related to a given signature.

    Returns
    -------
    V: array-like, shape=[n, n]
        A vertical vector.
    """
    return X * skew_diag_kernel


def flag_norm(H):
    """ Compute the norm of a tangent vector w.r.t. the flag canonical metric.

    Parameters
    ----------
    H: array-like, shape=[n, n]
        A horizontal vector.

    Returns
    -------
    ||H||_c: positive real
        The norm of H w.r.t the flag canonical metric.
    """
    return np.sqrt(1 / 2 * np.trace(H.T @ H))


def align_eucl(P, Q, signature_0):
    """ Use orthogonal Procrustes analysis as in [Szwagier2023, Theorem 2], to find the
    closest orthogonal matrix to P in the equivalence class of Q, in Frobenius distance.

    Parameters
    ----------
    P: array-like, shape=[n, n]
        An orthogonal matrix representing a flag.
    Q: array-like, shape=[n, n]
        An orthogonal matrix representing a flag.
    signature_0: tuple of ints
        The signature of the flag manifold.

    Returns
    -------
    Q_procrustes: array-like, shape=[n, n]
        An orthogonal matrix representing the closest orthogonal matrix to P in the
        equivalence class of Q, in Frobenius distance.
    """
    R_list = []
    angle_list = []
    for i in range(len(signature_0) - 1):
        P_i = P[:, signature_0[i]: signature_0[i + 1]]
        Q_i = Q[:, signature_0[i]: signature_0[i + 1]]
        U, s, Vh = svd(Q_i.T @ P_i)
        angle = np.arccos(np.clip(s, 0, 1))
        angle_list.append(angle)
        R_list.append(U @ Vh)
    R_procrustes = block_diag(*R_list)
    Q_procrustes = Q @ R_procrustes
    return Q_procrustes


def random_uniform_orbit(signature_0):
    r""" Sample a matrix from the uniform distribution on
    :math:`\mathcal{O}(\mathcal{I})`.

    Parameters
    ----------
    signature_0: tuple of ints
        The signature of the flag.

    Returns
    -------
    R: array-like, shape=[n, n]
        A random uniform matrix in :math:`\mathcal{O}(\mathcal{I})`.
    """
    flag_type = tuple(np.diff(signature_0))
    R_list = []
    for n_i in flag_type:
        if n_i == 1:
            R_list.append((2 * np.random.randint(2) - 1) * np.ones((1, 1)))
        else:
            R_list.append(ortho_group.rvs(dim=n_i))
    return block_diag(*R_list)


def random_uniform_special_orbit(signature_0):
    r""" Sample a matrix from the uniform distribution on
    :math:`\mathcal{SO}(\mathcal{I})`.

    Parameters
    ----------
    signature_0: tuple of ints
        The signature of the flag.

    Returns
    -------
    R: array-like, shape=[n, n]
        A random uniform matrix in :math:`\mathcal{SO}(\mathcal{I})`.
    """
    flag_type = tuple(np.diff(signature_0))
    R_list = []
    for n_i in flag_type:
        if n_i == 1:
            R_list.append(np.ones((1, 1)))
        else:
            R_list.append(special_ortho_group.rvs(dim=n_i))
    return block_diag(*R_list)


def switch_sign_cols(Q, signature_0):
    r""" Generate a list with all the equivalents of Q in the different connected
    components of :math:`\mathcal{O}(\mathcal{I})`, similarly as in [Ma2022].

    Parameters
    ----------
    Q: array-like, shape=[n, n]
        A flag.
    signature_0: tuple of ints
        The signature of the flag.

    Returns
    -------
    Q_list: list
        The list with all the equivalents of Q in the different connected components of
        :math:`\mathcal{O}(\mathcal{I})`.

    References
    ----------
    .. [Ma2022] X. Ma, M. Kirby, C.Peterson.
        “Self-organizing mappings on the flag manifold with applications to
        hyper-spectral image data analysis.” Neural Computing and Applications. 2022.
    """
    Q_list = []
    for j in range(len(signature_0[1:]) + 1):
        for cols_to_switch_sign in itertools.combinations(signature_0[:-1], j):
            Q_i = np.copy(Q)
            Q_i[:, cols_to_switch_sign] = - Q_i[:, cols_to_switch_sign]
            Q_list.append(Q_i)
    return Q_list


def flag_log(P, Q, signature_0, skew_diag_kernel, itermax=50, eps=1e-5):
    """ Algorithm for the approximation of the Riemannian logarithm on flag manifolds,
    using [Szwagier2023, Algorithm 1].

    Parameters
    ----------
    P: array-like, shape=[n, n]
        A flag.
    Q: array-like, shape=[n, n]
        A flag.
    signature_0: tuple of ints
        The signature of the flag.
    skew_diag_kernel: array-like, shape=[n, n]
        A binary mask with 1s on the diagonal blocks related to a given signature.
    itermax: int
        The maximal number of iterations.
    eps: positive real
        The endpoint error threshold under which we stop the algorithm.

    Returns
    -------
    H_: array-like, shape=[n, n]
        The horizontal vector approximating the Riemannian logarithm.
    flag_norm(H_): positive real
        Its norm, which is also the geodesic distance.
    err:
        The final endpoint error.
    count: int
        The number of iterations of the algorithm.
    time: positive real
        The total running time.
    err_list: list of real
        The history of endpoint errors along the iterations.
    """
    start = time()
    Q_ = np.copy(Q)
    err_list = []
    count = 0
    while True:
        X_ = logm(P.T @ Q_)
        H_ = proj_H(X_, skew_diag_kernel)
        M_ = P @ expm(H_)
        err = np.linalg.norm(
            projector_embedding(Q, signature_0) - projector_embedding(M_, signature_0),
            2)
        err_list.append(err)
        if (count >= itermax) or (err <= eps):
            break
        Q_ = align_eucl(M_, Q, signature_0)
        count += 1
    return H_, flag_norm(H_), err, count, time() - start, err_list


def flag_log_ma(P, Q, signature_0, skew_diag_kernel, itermax=50, eps=1e-5,
                try_all_orientations=True):
    """ Algorithm for the approximation of the Riemannian logarithm on flag manifolds,
    using [Ma2022, Algorithm 1].

    Parameters
    ----------
    P: array-like, shape=[n, n]
        A flag.
    Q: array-like, shape=[n, n]
        A flag.
    signature_0: tuple of ints
        The signature of the flag.
    skew_diag_kernel: array-like, shape=[n, n]
        A binary mask with 1s on the diagonal blocks related to a given signature.
    itermax: int
        The maximal number of iterations.
    eps: positive real
        The endpoint error threshold under which we stop the algorithm.
    try_all_orientations: bool
        Whether we try all the connected components in
        :math:`\mathcal{O}(\mathcal{I})`, like done in [Ma2022].

    Returns
    -------
    H_: array-like, shape=[n, n]
        The horizontal vector approximating the Riemannian logarithm.
    flag_norm(H_): positive real
        Its norm, which is also the geodesic distance.
    err:
        The final endpoint error.
    count: int
        The number of iterations of the algorithm.
    time: positive real
        The total running time.
    err_list: list of real
        The history of endpoint errors along the iterations.
    """
    start = time()
    Q_list = [Q]
    if try_all_orientations:
        Q_list = switch_sign_cols(Q, signature_0)
    H_opt, dist_opt, err_opt, count_opt = None, None, np.inf, None
    for Q in Q_list:
        try:
            G_ = random_uniform_special_orbit(signature_0)
            err_list = []
            count = 0
            while True:
                H_ = proj_H(logm((P.T @ Q) @ expm(-G_)), skew_diag_kernel)
                err = np.linalg.norm(
                    projector_embedding(Q, signature_0) - projector_embedding(
                        P @ expm(H_), signature_0), 2)
                err_list.append(err)
                if (count >= itermax) or (err <= eps):
                    break
                G_ = proj_V(logm(expm(-H_) @ (P.T @ Q)), skew_diag_kernel)
                count += 1
            if err < err_opt:
                H_opt, dist_opt, err_opt, count_opt, err_list_opt = H_, flag_norm(
                    H_), err, count, err_list
        except Exception:
            print("Algo exploded")
            continue
    return H_opt, dist_opt, err_opt, count_opt, time() - start, err_list_opt


signature = (1, 3, 5, 20, 100)
k, n = signature[-2:]
signature_0 = (0,) + signature
flag_type = tuple(np.diff(signature_0))
skew_diag_kernel = projection_kernel(signature_0)
path = "C:/Users/tszwagie/repos/Obsidian/Research Diary/Writing/GSI2023/figures/"

np.random.seed(42)
n_experiments = 10
itermax = 10
for n_exp, dist_true in enumerate([.2 * np.pi, .5 * np.pi, 1 * np.pi]):
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    results_mean = np.zeros((2, 5))
    results_std = np.zeros((2, 5))
    results_align = np.zeros((n_experiments, 5))
    results_proj = np.zeros((n_experiments, 5))
    for i in range(n_experiments):
        print(i)
        P = ortho_group.rvs(dim=n)
        H_true = proj_H(SkewSymmetricMatrices(n).random_point(), skew_diag_kernel)
        H_true = H_true / flag_norm(H_true) * dist_true
        Q_aligned = P @ expm(H_true)
        R_true = random_uniform_orbit(signature_0)
        Q = Q_aligned @ R_true

        Q_procrustes = align_eucl(P, Q, signature_0)
        H, dist, err, count, dt, err_list = flag_log(P, Q_procrustes, signature_0,
                                                     skew_diag_kernel, itermax=itermax,
                                                     eps=1e-5)
        results_align[i] = np.array(
            [err, abs(dist - dist_true), np.linalg.norm(H - H_true, 2), count, dt])
        ax.plot(np.arange(0, count + 1), err_list, color='tab:red',
                label='Alignment' if i == 0 else None)

        H, dist, err, count, dt, err_list = flag_log_ma(P, Q, signature_0,
                                                        skew_diag_kernel,
                                                        itermax=itermax, eps=1e-5,
                                                        try_all_orientations=True)
        results_proj[i] = np.array(
            [err, abs(dist - dist_true), np.linalg.norm(H - H_true, 2), count, dt])
        ax.plot(np.arange(0, count + 1), err_list, color='tab:blue', ls='--',
                label='Alternate projections' if i == 0 else None)

    ax.set_yscale('log')
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Endpoint error')
    plt.legend()
    plt.savefig(path + f"Exp{n_exp+1}_(1,3,5,20,100).png", dpi='figure', format='png', transparent=True)
    plt.savefig(path + f"Exp{n_exp+1}_(1,3,5,20,100).pdf", dpi='figure', format='pdf', transparent=True)
    plt.savefig(path + f"Exp{n_exp+1}_(1,3,5,20,100).svg", dpi='figure', format='svg', transparent=True)

    results_mean[0] = np.mean(results_align, axis=0)
    results_mean[1] = np.mean(results_proj, axis=0)
    results_std[0] = np.std(results_align, axis=0)
    results_std[1] = np.std(results_proj, axis=0)
    df_results_mean = pd.DataFrame(results_mean,
                                   index=['Alignment', 'Alternate projections'],
                                   columns=['endpoint error', 'dist error',
                                            'tangent error', 'iterations', 'time'])
    df_results_std = pd.DataFrame(results_std,
                                  index=['Alignment', 'Alternate projections'],
                                  columns=['endpoint error', 'dist error',
                                           'tangent error', 'iterations', 'time'])
    print(df_results_mean.to_latex(float_format="%.1e"))
    print(df_results_std.to_latex(float_format="%.1e"))
