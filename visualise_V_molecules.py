import matplotlib.pyplot as plt
from causalset import CausalSet
from causalsetfunctions import compute_spacetimecuts_uniform_Rindler,  compute_spacetimecuts_tube


def main():
    # boundsArray, adjusted_rho, l = compute_spacetimecuts_uniform_Rindler(
    #     d=3, rho2=1000, N_max=8000, b=3)
    # C = CausalSet(sprinkling_density=adjusted_rho, dimension=3,
    #               BHtype='Rindler', bounds=boundsArray)

    boundsArray, rho = compute_spacetimecuts_tube(
        d=3, rho2=10, N_max=1000, b=3)
    C = CausalSet(sprinkling_density=rho, dimension=3,
                  BHtype='Dynamic', sprinkling='Tube', T=1, bounds=boundsArray)

    return C, boundsArray


def main2():

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
    C, boundsArray = main()
    main2(C, boundsArray)
