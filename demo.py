import numpy as np
import matplotlib.pyplot as plt
from causalset import CausalSet
plt.rc('font', family='Arial', size=16)

def main():
    C = CausalSet(sprinkling_density=10, dimension=3, BHtype='Rindler', bounds=np.array([[-0.5, 0.5] for i in range(3)]))
    C.find_linkmatrix()
    coordinates = np.array([x.coordinates for x in C.ElementList])

    # plt.scatter(coordinates[:, 1], coordinates[:, 0])
    # print(C.LinkMatrix)
    # for i in range(len(C.LinkMatrix)):
    #     for j in range(len(C.LinkMatrix[i])):
    #         if C.LinkMatrix[i][j] == 1:
    #             plt.plot([coordinates[i][1], coordinates[j][1]], [coordinates[i][0], coordinates[j][0]], 'r', zorder=0)

    # 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.view_init(15, 45)
    ax.scatter3D(coordinates[:, 0], coordinates[:, 1], coordinates[:, 2], 'o', color='black', s=100)

    # print(C.LinkMatrix)
    for i in range(len(C.LinkMatrix)):
        for j in range(len(C.LinkMatrix[i])):
            if C.LinkMatrix[i][j] == 1:
                ax.plot([coordinates[i][0], coordinates[j][0]], [coordinates[i][1], coordinates[j][1]], [coordinates[i][2], coordinates[j][2]], color='#002060')

    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_zlabel('$z$')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()