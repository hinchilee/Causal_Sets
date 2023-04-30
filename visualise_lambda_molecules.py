import matplotlib.pyplot as plt
from causalset import CausalSet
from causalsetfunctions import compute_spacetimecuts_uniform_Rindler, compute_spacetimecuts_tube
import numpy as np


def main():
    np.random.seed(10)
    # boundsArray, adjusted_rho, l = compute_spacetimecuts_uniform_Rindler(
    #     d=4, rho2=1000, N_max=10000, b=3)
    # C = CausalSet(sprinkling_density=adjusted_rho, dimension=4,
    #               BHtype='Rindler', bounds=boundsArray)

    boundsArray, rho = compute_spacetimecuts_tube(
        d=4, rho2=10, N_max=18000, b=1.9)
    C = CausalSet(sprinkling_density=rho, dimension=4,
                  BHtype='Dynamic', sprinkling='Tube', T=1, bounds=boundsArray)
    return C, boundsArray


def main2(C, boundsArray):

    C.find_molecules()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # ax.axes.set_xlim3d(left=boundsArray[0][0], right=boundsArray[0][1])
    # ax.axes.set_ylim3d(bottom=boundsArray[1][0], top=boundsArray[1][1])
    # ax.axes.set_zlim3d(bottom=boundsArray[2][0], top=boundsArray[2][1])

    # ax.axes.set_xlim3d(left=boundsArray[2], right=boundsArray[3])
    # ax.axes.set_ylim3d(bottom=-boundsArray[1], top=boundsArray[1])
    # ax.axes.set_zlim3d(bottom=-boundsArray[1], top=boundsArray[1])

    ax.axes.set_xlim3d(left=-1, right=1)
    ax.axes.set_ylim3d(bottom=-1, top=1)
    ax.axes.set_zlim3d(bottom=-1, top=1)

    # C.find_molecules()
    # print(C.MoleculesList)
    # print(C.VElementsLabelsList)

    t = []
    x = []
    y = []
    for i in C.MoleculesList:
        t.append(C.ElementList[i].coordinates[1])
        x.append(C.ElementList[i].coordinates[2])
        y.append(C.ElementList[i].coordinates[3])

    for i in range(len(C.MoleculesList) // 2):
        ax.plot([t[2*i], t[2*i+1]], [x[2*i], x[2*i+1]],
                [y[2*i], y[2*i+1]], color='black')

    ax.scatter([t[i] for i in range(len(t)) if i % 2 == 0],
               [x[i] for i in range(len(x)) if i % 2 == 0],
               [y[i] for i in range(len(y)) if i % 2 == 0], s=80, alpha=1, label='Maximal-but-ones')

    ax.scatter([t[i] for i in range(len(t)) if i % 2 == 1],
               [x[i] for i in range(len(x)) if i % 2 == 1],
               [y[i] for i in range(len(y)) if i % 2 == 1], s=80, alpha=1, label='Maximals')

    params = {'font.family': 'Times New Roman',
              'xtick.labelsize': 27,
              'ytick.labelsize': 27,
              'axes.labelsize': 40,
              'font.size': 30,
              'axes.labelpad': 1
              }
    plt.rcParams.update(params)
    ax.set_xlabel('x', fontsize=30)
    ax.set_ylabel('y', fontsize=30)
    ax.set_zlabel('z', fontsize=30)
    ax.set_xticks([-1, -0.5, 0, 0.5, 1])
    ax.set_yticks([-1, -0.5, 0, 0.5, 1])
    ax.set_zticks([-1, -0.5, 0, 0.5, 1])
    plt.legend()
    plt.show()


def main3(C, boundsArray):

    V, b2 = C.find_Vmolecules()
    # print(C.VElementsLabelsList)

    t = []
    x = []
    y = []
    for i in C.VElementsLabelsList:
        t.append(C.ElementList[i].coordinates[1])
        x.append(C.ElementList[i].coordinates[2])
        y.append(C.ElementList[i].coordinates[3])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.axes.set_xlim3d(left=boundsArray[1][0], right=boundsArray[1][1])
    ax.axes.set_ylim3d(bottom=boundsArray[2][0], top=boundsArray[2][1])
    ax.axes.set_zlim3d(bottom=boundsArray[3][0], top=boundsArray[3][1])

    for i in range(len(C.VElementsLabelsList) // 3):
        ax.plot([t[3*i], t[3*i+2]], [x[3*i], x[3*i+2]],
                [y[3*i], y[3*i+2]], color='black')
        ax.plot([t[3*i+1], t[3*i+2]], [x[3*i+1], x[3*i+2]],
                [y[3*i+1], y[3*i+2]], color='black')

    ax.scatter(t, x, y)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()


if __name__ == "__main__":
    #C, boundsArray = main()
    # lambda-molecules
    #main2(C, boundsArray)
    # v moleculues
    main3(C, boundsArray)
