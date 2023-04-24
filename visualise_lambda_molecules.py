import matplotlib.pyplot as plt
from causalset import CausalSet
from causalsetfunctions import compute_spacetimecuts_uniform_Rindler, compute_spacetimecuts_tube

def main():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # boundsArray, adjusted_rho, l, adjusted_l = compute_spacetimecuts_uniform_Rindler(d=3, rho=10000, N_max=8000, b=3)
    # C = CausalSet(sprinkling_density=adjusted_rho, dimension=3, BHtype='Rindler', bounds = boundsArray)

    # ax.axes.set_xlim3d(left=boundsArray[0][0], right=boundsArray[0][1])
    # ax.axes.set_ylim3d(bottom=boundsArray[1][0], top=boundsArray[1][1])
    # ax.axes.set_zlim3d(bottom=boundsArray[2][0], top=boundsArray[2][1])

    boundsArray, rho = compute_spacetimecuts_tube(d=3, rho2=10, N_max=8000, b=3)
    C = CausalSet(sprinkling_density=rho, dimension=3, BHtype='Dynamic', sprinkling='Tube', T=1, bounds=boundsArray)

    ax.axes.set_xlim3d(left=boundsArray[2], right=boundsArray[3])
    ax.axes.set_ylim3d(bottom=-boundsArray[1], top=boundsArray[1])
    ax.axes.set_zlim3d(bottom=-boundsArray[1], top=boundsArray[1])


    C.find_molecules()
    # print(C.MoleculesList)
    # print(C.VElementsLabelsList)

    t = []
    x = []
    y = []
    for i in C.MoleculesList:
        t.append(C.ElementList[i].coordinates[0])
        x.append(C.ElementList[i].coordinates[1])
        y.append(C.ElementList[i].coordinates[2])

    

    for i in range(len(C.MoleculesList) // 2):
        ax.plot([t[2*i], t[2*i+1]], [x[2*i], x[2*i+1]], [y[2*i], y[2*i+1]], color='black')

    ax.scatter([t[i] for i in range(len(t)) if i % 2 == 0],
                [x[i] for i in range(len(x)) if i % 2 == 0],
                [y[i] for i in range(len(y)) if i % 2 == 0])

    ax.scatter([t[i] for i in range(len(t)) if i % 2 == 1],
                [x[i] for i in range(len(x)) if i % 2 == 1],
                [y[i] for i in range(len(y)) if i % 2 == 1])

    ax.set_xlabel('t')
    ax.set_ylabel('x')
    ax.set_zlabel('y')
    plt.show()

if __name__ == "__main__":
    main()