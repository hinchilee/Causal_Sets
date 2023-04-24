import matplotlib.pyplot as plt
from causalset import CausalSet
from causalsetfunctions import compute_spacetimecuts_uniform_Rindler

def main():
    boundsArray, adjusted_rho, l, adjusted_l = compute_spacetimecuts_uniform_Rindler(d=4, rho=30000, N_max=10000, b=3)
    C = CausalSet(sprinkling_density=adjusted_rho, dimension=4, BHtype='Rindler', bounds = boundsArray)
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
        ax.plot([t[3*i], t[3*i+2]], [x[3*i], x[3*i+2]], [y[3*i], y[3*i+2]], color='black')
        ax.plot([t[3*i+1], t[3*i+2]], [x[3*i+1], x[3*i+2]], [y[3*i+1], y[3*i+2]], color='black')

    ax.scatter(t, x, y)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()

if __name__ == "__main__":
    main()