import numpy as np
import matplotlib.pyplot as plt

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

def approximates(a, b):
    err = 10**(-6)
    return b - err < a and a < b + err

def Haldane_zigzagx(Nx, Ny):
    # Nx, Ny: even, 六角格子を各方向に並べる数
    sites = []
    xs = np.arange(0, 3 * Nx / 2 + 1, 1 / 2)

    ys_0_2_4 = np.arange(np.sqrt(3) / 2, (Ny + 3 / 2) * np.sqrt(3), np.sqrt(3))
    ys_1_3_5 = np.arange(0, (Ny + 1) * np.sqrt(3), np.sqrt(3))

    for x in xs:
        if x % 3 == 0:
            for y in ys_0_2_4:
                sites += [[x, y, "a"]]
        elif (x - 1 / 2) % 3 == 0:
            for y in ys_1_3_5:
                sites += [[x, y, "b"]]
        elif (x - 1) % 3 == 0:
            for y in ys_0_2_4:
                sites += [[x, y, "c"]]
        elif (x - 3 / 2) % 3 == 0:
            for y in ys_1_3_5:
                sites += [[x, y, "a"]]
        elif (x - 2) % 3 == 0:
            for y in ys_0_2_4:
                sites += [[x, y, "b"]]
        elif (x - 5 / 2) % 3 == 0:
            for y in ys_1_3_5:
                sites += [[x, y, "c"]]

    return sites

# ホッピングインデックスの獲得
def get_hopping_indices(sites):
    err = 10**(-6)
    nearests = []
    triangles = []
    couplings = []
    next_nearests = []

    for i in range(len(sites)):
        for j in range(i, len(sites)):
            xi, yi, si = sites[i][0], sites[i][1], sites[i][2]
            xj, yj, sj = sites[j][0], sites[j][1], sites[j][2]
            xij = xi - xj
            yij = yi - yj

            distance = np.sqrt(xij**2 + yij**2)
            if approximates(distance, np.sqrt(3)) and (si == "c" and sj == "c"):
                xij = xij / np.sqrt(3)
                yij = yij / np.sqrt(3)
                triangles.append([i, j, xij, yij])
            elif approximates(distance, 1) and (si == "c" or sj == "c"):
                couplings.append([i, j, xij, yij])
            elif approximates(distance, 1):
                nearests.append([i, j, xij, yij])
            elif approximates(distance, np.sqrt(3)):
                xij = xij / np.sqrt(3)
                yij = yij / np.sqrt(3)
                next_nearests.append([i, j, xij, yij])

    return nearests, triangles, couplings, next_nearests

def make_unitvector_a_rij():
    v0 = [0, -1]
    v1 = [np.sqrt(3) / 2, -1 / 2]
    v2 = [np.sqrt(3) / 2, 1 / 2]
    v3 = [1, 0]
    v4 = [-np.sqrt(3), 1 / 2]
    v5 = [-np.sqrt(3) / 2, -1 / 2]

    return [v0, v1, v2, v3, v4, v5]

def identify_direction(vs, sublattice, next_nearest):
    xij = next_nearest[2]
    yij = next_nearest[3]

    for i in range(len(vs)):
        if (sublattice == "a") and (i == 0 or i == 2 or i == 4) and (approximates(xij, vs[i][0]) and approximates(yij, vs[i][1])):
            return "+"
        elif (sublattice == "a") and (i == 1 or i == 3 or i == 5) and (approximates(xij, vs[i][0]) and approximates(yij, vs[i][1])):
            return "-"
        elif (sublattice == "b") and (i == 1 or i == 3 or i == 5) and (approximates(xij, vs[i][0]) and approximates(yij, vs[i][1])):
            return "+"
        else:
            return "-"

def add_phase_info(sites, next_nearests):
    vs = make_unitvector_a_rij()
    sublattices = [sites[i][2] for i in range(len(sites))]

    for next_nearest in next_nearests:
        i = next_nearest[0]
        sign = identify_direction(vs, sublattices[i], next_nearest)
        next_nearest += [sign]

    return next_nearests

def make_plot_hop(sites, nearests, triangles, couplings, next_nearests):
    fig, ax = plt.subplots(figsize=(10,10))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    fig.patch.set_alpha(0)
    ax.patch.set_alpha(0)

    for X in sites:
        colors_dict = {"a": "orange", "b": "green", "c": "navy"}
        marker_dict = {"a": "o", "b": "o", "c": "^"}
        ax.scatter(X[0], X[1], zorder=3, marker = marker_dict[X[2]], c=colors_dict[X[2]])

    for hop in nearests:
        xi = sites[hop[0]][0]
        yi = sites[hop[0]][1]
        xitoj = -hop[2]
        yitoj = -hop[3]
        ax.plot([xi, xi + xitoj], [yi, yi + yitoj], ls = "-", c = "black")

    for hop in triangles:
        xi = sites[hop[0]][0]
        yi = sites[hop[0]][1]
        xitoj = -hop[2] * np.sqrt(3)
        yitoj = -hop[3] * np.sqrt(3)
        ax.plot([xi, xi + xitoj], [yi, yi + yitoj], ls = "--", c = "aqua")

    for hop in couplings:
        xi = sites[hop[0]][0]
        yi = sites[hop[0]][1]
        xitoj = -hop[2]
        yitoj = -hop[3]
        ax.plot([xi, xi + xitoj], [yi, yi + yitoj], ls = ":", c = "gold")

    for hop in next_nearests:
        xi = sites[hop[0]][0]
        yi = sites[hop[0]][1]
        xitoj = -hop[2] * np.sqrt(3)
        yitoj = -hop[3] * np.sqrt(3)
        ax.plot([xi, xi + xitoj], [yi, yi + yitoj], ls = "-.", c = "gray")

    plt.show()

def make_dense_Hamiltonian(sites, nearests, triangles, couplings, next_nearests):
    N = len(sites)
    H = np.zeros((N, N), dtype="complex128")

    # on-site potential
    m_dict = {"a": -m, "b": m, "c": m2}
    for k in range(N):
        i = k
        j = k
        Hij = m_dict[sites[k][2]]

        H[i, j] = Hij

    # nearest hopping
    for hop in nearests:
        i = hop[0]
        j = hop[1]
        Hij = -t1

        H[i, j] = Hij
        H[j, i] = Hij

    # triangle hopping
    for hop in triangles:
        i = hop[0]
        j = hop[1]
        Hij = -t2

        H[i, j] = Hij
        H[j, i] = Hij

    # hexagonal-triangle hopping
    for hop in couplings:
        i = hop[0]
        j = hop[1]
        Hij = -t3

        H[i, j] = Hij
        H[j, i] = Hij

    # next-nearest hopping
    phi_dict = {"+": phi, "-": -phi}
    for hop in next_nearests:
        i = hop[0]
        j = hop[1]
        Hij = -tc * np.exp(1j * phi_dict[hop[4]])

        H[i, j] = Hij
        H[j, i] = Hij.conjugate()

    return H

def make_spectral_localizer(x, y, kappa, X, Y, H):
    # E = 0で固定
    Ntot = len(X)
    X = np.diag(X - x)
    Y = np.diag(Y - y)
    L = np.zeros((2 * Ntot, 2 * Ntot), dtype="complex128")

    L[:Ntot, :Ntot] = H
    L[Ntot:, :Ntot] = kappa * (X + 1j * Y)
    L[:Ntot, Ntot:] = kappa * (X - 1j * Y)
    L[Ntot:, Ntot:] = -H

    return L

def local_Chern_number(x, y, X, Y, H):
    #_, D = LDLT_decomposition(L) こっちはゼロ除算が出てしまって上手くいかない（対角化の結果と変わってしまう）
    L = make_spectral_localizer(x, y, kappa, X, Y, H)
    D, _ = np.linalg.eigh(L)

    CL = 0
    gap = 10
    for d in D:
        if d > 0:
            CL += 1
        elif d < 0:
            CL -= 1

        if abs(d) < abs(gap):
            gap = d

    CL /= 2

    return CL, gap

def gap_of_x(x, y, ind_y, gaps):
    gaps = gaps.reshape(len(xs), len(ys))

    return gaps[:, ind_y]

def make_plot_CL(x, y, data):
    # data = array([CLi])
    X, Y = np.meshgrid(x, y)
    Z = data.reshape((len(x), len(y))).T

    fig, ax = plt.subplots(figsize=(13,10))
    im = ax.imshow(Z, interpolation='gaussian', cmap="plasma", origin='lower', vmax=Z.max(), vmin=Z.min())
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")

    plt.colorbar(im, label='Local Chern number')

    plt.show()

def make_plot_gap(x, y, data):
    # data = array([gapi])
    X, Y = np.meshgrid(x, y)
    Z = np.abs(data.reshape((len(x), len(y))).T)

    fig, ax = plt.subplots(figsize=(13,10))
    im = ax.imshow(Z, interpolation='gaussian', cmap="plasma", origin='lower', vmax=1, vmin=0)
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")

    plt.colorbar(im, label='Localizer gap')

    plt.show()

def make_plot_spectrum(x, spectrum):
    fig, ax = plt.subplots(figsize = (10, 10))
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"Energy $E$ / $t_1$")
    ax.set_xlim(-1, 18.686533479473212)
    ax.set_ylim(-0.3, 0.3)
    ax.minorticks_on()

    ax.plot(x, spectrum, c="gray")
    ax.hlines([0.], xmin=-2, xmax=19, lw=1, colors="black")

    plt.show()



# 実行
Nx = 10
Ny = 10
sites = Haldane_zigzagx(Nx, Ny)
nearests, triangles, couplings, next_nearests = get_hopping_indices(sites)
next_nearests = add_phase_info(sites, next_nearests)
H = make_dense_Hamiltonian(sites, nearests, triangles, couplings, next_nearests)
X = np.array([sites[i][0] for i in range(len(sites))])
Y = np.array([sites[i][1] for i in range(len(sites))])

kappa = 1
xs = np.arange(-1, 16.5, 0.5)
ys = np.arange(-1, 18.686533479473212, 0.5)

data_CH = []
data_gap = []
cnt = 0
for x in xs:
    for y in ys:
        CL, gap = local_Chern_number(x, y, X, Y, H)
        data_CH += [CL]
        data_gap += [gap]
        cnt += 1
        print(cnt)

data_CH = np.array(data_CH)
data_gap = np.array(data_gap)

make_plot_CL(xs, ys, data_CH)
make_plot_gap(xs, ys, data_gap)

spectrum = gap_of_x(x, y, 14, data_gap)
make_plot_spectrum(xs, spectrum)
