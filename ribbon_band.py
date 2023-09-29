import numpy as np

N = 5
Nx = 2 * N
ky = np.linspace(0, 2 * np.pi / np.sqrt(3), 100)
# on-site potential at sublattice a(-), b(+)
m = 0
# on-site potential at triangle lattice
m2 = -0.35
# nearest hopping
t1 = 1
# triangle lattice hopping
t2 = 0.2
# hexagonal-triangle hopping
t3 = 0.3
# next nearest hopping
tc = 0.5
phi = np.pi / 2

def make_xinds():
    x_inds = {"a0": [i for i in range(N + 1)],
              "a1": [i for i in range(N + 1, 2 * N + 1)],
              "b0": [i for i in range(2 * N + 1, 3 * N + 2)],
              "b1": [i for i in range(3 * N + 2, 4 * N + 2)],
              "c0": [i for i in range(4 * N + 2, 5 * N + 2)],
              "c1": [i for i in range(5 * N + 2, 6 * N + 2)]}

    return x_inds

def make_Hamiltonian_xky(x_inds, ky):
    H = np.zeros((6 * N + 2, 6 * N + 2))

    # on-site potential
    for k in range(6 * N + 2):
        i = k
        j = k
        if k < 2 * N + 1:
            Hij = -m
        elif k < 4 * N + 2:
            Hij = m
        else:
            Hij = m2

        H[i, j] += Hij

    # nearest hopping
    for k in range(N + 1):
        i = x_inds["a0"][k]
        j = x_inds["b0"][k]
        Hij = -2 * t1 * np.cos(np.sqrt(3) * ky / 2)

        H[i, j] += Hij
        H[j, i] += Hij

        if k == N:
            break

        i = x_inds["a1"][k]
        j = x_inds["b1"][k]
        Hij = -2 * t1 * np.cos(np.sqrt(3) * ky / 2)

        H[i, j] += Hij
        H[j, i] += Hij

        i = x_inds["b0"][k]
        j = x_inds["a1"][k]
        Hij = -t1

        H[i, j] += Hij
        H[j, i] += Hij

        i = x_inds["b1"][k]
        j = x_inds["a0"][k + 1]
        Hij = -t1

        H[i, j] += Hij
        H[j, i] += Hij

    # triangle hopping
    for k in range(N):
        i = x_inds["c0"][k]
        j = x_inds["c0"][k]
        Hij = -2 * t2 * np.cos(np.sqrt(3) * ky)

        H[i, j] += Hij

        i = x_inds["c1"][k]
        j = x_inds["c1"][k]

        H[i, j] += Hij

        i = x_inds["c0"][k]
        j = x_inds["c1"][k]
        Hij = -2 * t2 * np.cos(np.sqrt(3) * ky / 2)

        H[i, j] += Hij
        H[j, i] += Hij

        if k == N - 1:
            break

        i = x_inds["c0"][k + 1]

        H[i, j] += Hij
        H[j, i] += Hij

    # honeycomb-triangle hopping
    for k in range(N):
        ijs = [(x_inds["c0"][k], x_inds["b1"][k]),
               (x_inds["c0"][k], x_inds["a0"][k]),
               (x_inds["c1"][k], x_inds["b0"][k + 1]),
               (x_inds["c1"][k], x_inds["a1"][k])]
        Hij = -t3
        for ij in ijs:
            i, j = ij[0], ij[1]
            H[i, j] += Hij
            H[j, i] += Hij

        ijs = [(x_inds["c0"][k], x_inds["b0"][k]),
               (x_inds["c0"][k], x_inds["a1"][k]),
               (x_inds["c1"][k], x_inds["b1"][k]),
               (x_inds["c1"][k], x_inds["a0"][k + 1])]
        Hij = -2 * t3 * np.cos(np.sqrt(3) * ky / 2)

        for ij in ijs:
            i, j = ij[0], ij[1]
            H[i, j] += Hij
            H[j, i] += Hij

    # next nearest hopping
    for k in range(N):
        ij_Hijs = [(x_inds["a0"][k], x_inds["a0"][k], -2 * tc * np.cos(np.sqrt(3) * ky + phi)),
                   (x_inds["a1"][k], x_inds["a1"][k], -2 * tc * np.cos(np.sqrt(3) * ky + phi)),
                   (x_inds["a0"][k], x_inds["a1"][k], -2 * tc * np.cos(np.sqrt(3) * ky / 2 - phi)),
                   (x_inds["a1"][k], x_inds["a0"][k + 1], -2 * tc * np.cos(np.sqrt(3) * ky / 2 - phi)),
                   (x_inds["b0"][k], x_inds["b0"][k], -2 * tc * np.cos(np.sqrt(3) * ky - phi)),
                   (x_inds["b1"][k], x_inds["b1"][k], -2 * tc * np.cos(np.sqrt(3) * ky - phi)),
                   (x_inds["b0"][k], x_inds["b1"][k], -2 * tc * np.cos(np.sqrt(3) * ky / 2 + phi)),
                   (x_inds["b1"][k], x_inds["b0"][k + 1], -2 * tc * np.cos(np.sqrt(3) * ky / 2 + phi))]

        for ij_Hij in ij_Hijs:
            i, j, Hij = ij_Hij[0], ij_Hij[1], ij_Hij[2]

            H[i, j] += Hij
            if i != j:
                H[j, i] += Hij

    # 右端のオンサイトが足されていないので、足す。
    ij_Hijs = [(x_inds["a0"][N], x_inds["a0"][N], -2 * tc * np.cos(np.sqrt(3) * ky + phi)),
               (x_inds["b0"][N], x_inds["b0"][N], -2 * tc * np.cos(np.sqrt(3) * ky - phi))]

    for ij_Hij in ij_Hijs:
        i, j, Hij = ij_Hij[0], ij_Hij[1], ij_Hij[2]

        H[i, j] += Hij

    return H

def Energy(kys):
    Es = []
    x_inds = make_xinds()

    for ky in kys:
        H = make_Hamiltonian_xky(x_inds, ky)
        E, _ = np.linalg.eigh(H)
        Es += [E]

    return Es

import matplotlib.pyplot as plt
plt.rcParams["font.size"] = 24

def make_plot_energy(kys, Es):
    kys = np.sqrt(3) * kys / np.pi
    fig, ax = plt.subplots(figsize = (10, 10))
    ax.set_xlabel(r"Wave number $\sqrt{3} k_y$ / $\pi$")
    ax.set_ylabel(r"Energy $E$ / $t_1$")
    ax.set_xlim(0, 2)
    ax.set_ylim(-4, 4)
    ax.minorticks_on()

    for i in range(len(kys)):
        k = kys[i]
        for E in Es[i]:
            ax.scatter(k, E, zorder=3, c="blue")

    plt.show()



# 実行
Es = Energy(ky)
make_plot_energy(ky, Es)
