import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import jaccard_score
from causalset import CausalSet

C = CausalSet(dimension=3, number_of_points=10)
coordinates = np.array([x.coordinates for x in C.List])

plt.scatter(coordinates[:, 1], coordinates[:, 0])
print(C.LinkMatrix)
for i in range(len(C.LinkMatrix)):
    for j in range(len(C.LinkMatrix[i])):
        if C.LinkMatrix[i][j] == 1:
            plt.plot([coordinates[i][1], coordinates[j][1]], [coordinates[i][0], coordinates[j][0]], 'r')

# 3D plot
# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# ax.scatter3D(coordinates[:, 0], coordinates[:, 1], coordinates[:, 2], 'o')

# print(C.LinkMatrix)
# for i in range(len(C.LinkMatrix)):
#     for j in range(len(C.LinkMatrix[i])):
#         if C.LinkMatrix[i][j] == 1:
#             ax.plot([coordinates[i][0], coordinates[j][0]], [coordinates[i][1], coordinates[j][1]], [coordinates[i][2], coordinates[j][2]], 'r')

plt.show()