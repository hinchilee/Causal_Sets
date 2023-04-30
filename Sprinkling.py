
import numpy as np


def Sprinkling_Uniform(dimension=2, number_of_points=10, bounds=np.array([[0, 1], [0, 1]])):
    # Generate a list of n points of (t, x1, x2, ..., x_d) randomly in the bounds
    # bounds is a nparray of (low, high) bounds, first ele is time, the remaining is space

    if len(bounds) != dimension:
        raise ValueError(
            'The dimension of bounds does not match the bounds of the sprinkling dimension!')

    sprinkledCoords: np.array = np.random.rand(number_of_points, dimension)

    # Scaling the distribution to fit the time bounds
    sprinkledCoords[:, 0] *= (bounds[0][1] - bounds[0][0])
    sprinkledCoords[:, 0] += (bounds[0][0])

    # Scaling the distribution to fit the space bounds
    for xn in range(1, dimension):
        sprinkledCoords[:, xn] *= (bounds[xn][1] - bounds[xn][0])
        sprinkledCoords[:, xn] += (bounds[xn][0])

    return sprinkledCoords


def Sprinkling_Bicone(dimension=2, number_of_points=10):

    Time = list()
    Space = list()
    T = 1
    R_0 = 1
    for i in range(number_of_points):
        x = np.random.uniform(0, 1)
        sign = np.random.uniform(0, 1)
        t = T*(1-x**(1/dimension))
        if sign > 0.5:
            Time.append(t)
        else:
            Time.append(-t)
        R_t = R_0 * x ** (1/dimension)
        v_unnormalised = np.random.normal(0, 1, size=dimension - 1)

        def normalise(v):
            norm = np.linalg.norm(v)
            if norm == 0:
                return v
            return v/norm

        v = normalise(v_unnormalised)
        k = np.random.uniform(0, 1)
        v_f = v*R_t*k**(1/(dimension - 1))
        Space.append(v_f)

    Time = np.array(Time)
    Space = np.array(Space)
    Coords = np.concatenate((Time.reshape(number_of_points, 1), Space), axis=1)

    # Manually add top and bottom points
    top_point = np.concatenate(
        ([1], np.zeros(shape=(dimension-1, ), dtype=int)))
    bottom_point = np.concatenate(
        ([-1], np.zeros(shape=(dimension-1, ), dtype=int)))
    Coords = np.concatenate((Coords, top_point.reshape(1, dimension)), axis=0)
    Coords = np.concatenate(
        (Coords, bottom_point.reshape(1, dimension)), axis=0)

    return Coords


def Sprinkling_Tube(dimension=2, number_of_points=10, bounds=[0.3, 1, 0, 1]):

    R_min, R_max, T_min, T_max = bounds
    Time = list()
    Space = list()
    n = dimension - 1
    x_min = (R_min / R_max)**(n)
    for i in range(number_of_points):
        x = np.random.uniform(0, 1)
        # uniform distribution between k_min^n to 1
        xprime = (1-x_min)*x + x_min
        R_t = R_max * xprime ** (1/(n))
        v_unnormalised = np.random.normal(0, 1, size=dimension - 1)

        def normalise(v):
            norm = np.linalg.norm(v)
            if norm == 0:
                return v
            return v/norm
        v = normalise(v_unnormalised)
        v_f = v*R_t
        Space.append(v_f)

        t2 = np.random.uniform(0, 1)
        # uniform distribution between T_min and T_max
        t = (T_max-T_min)*t2 + T_min
        Time.append(t)

    Time = np.array(Time)
    Space = np.array(Space)
    Coords = np.concatenate((Time.reshape(number_of_points, 1), Space), axis=1)
    return Coords


if __name__ == "__main__":
    # coords = Sprinkling_Bicone(dimension = 4, number_of_points = 10)
    np.random.seed(20)

    from scipy.optimize import curve_fit
    import matplotlib.pyplot as plt

    import matplotlib
    matplotlib.rcParams.update(matplotlib.rcParamsDefault)

    params = {'font.size': 25,
              'font.family': 'Times New Roman',
              'axes.labelsize': 23,
              'legend.fontsize': 15,
              'xtick.labelsize': 19,
              'ytick.labelsize': 19,
              'figure.figsize': [10.5, 6.5],
              'axes.prop_cycle': plt.cycler(color=plt.rcParams['axes.prop_cycle'].by_key()['color']
                                            + ['magenta'])
              }
    plt.rcParams.update(params)

    xLinspace = np.linspace(0, 1, 100)
    plt.figure(figsize=(6, 5))
    a = Sprinkling_Tube(dimension=4,
                        number_of_points=30000,
                        bounds=[0.2, 1, 0, 1])

    r = list()
    for row in range(a.shape[0]):
        r.append(np.linalg.norm(a[row][1:]))

        # bins = np.linspace(0, 1, 21)

    freq, binsPositions = np.histogram(r, bins=21)
    binMiddles = list()
    for i in range(len(binsPositions)-1):
        binMiddles.append((binsPositions[i] + binsPositions[i+1])/2)
    plt.hist(r, bins=binsPositions, label='Sprinkled Points 3+1',
             alpha=0.7, color='green')
    binMiddles = np.array(binMiddles)

    def funcx3(x, a):
        return a*x**2

    popt, pcov = curve_fit(funcx3, binMiddles, freq)
    plt.plot(xLinspace, funcx3(xLinspace, *popt),
             label=r'$r^2$ fit', color='orange')

    b = Sprinkling_Tube(dimension=3,
                        number_of_points=30000,
                        bounds=[0.2, 1, 0, 1])

    r = list()
    for row in range(b.shape[0]):
        r.append(np.linalg.norm(b[row][1:]))

    freq, binsPositions = np.histogram(r, bins=21)
    binMiddles = list()
    for i in range(len(binsPositions)-1):
        binMiddles.append((binsPositions[i] + binsPositions[i+1])/2)
    plt.hist(r, bins=binsPositions, label='Sprinkled Points 2+1',
             alpha=0.4, color='blue')
    binMiddles = np.array(binMiddles)

    def funcx2(x, a):
        return a*x**1

    popt, pcov = curve_fit(funcx2, binMiddles, freq)
    plt.yticks([1000, 2000, 3000])
    plt.xticks(np.arange(0, 1.2, 0.2))
    plt.plot(xLinspace, funcx2(xLinspace, *popt),
             label=r'$r$ fit', color='red', alpha=0.8)

    plt.xlim(0, 1)
    plt.xlabel(r'Radius from central axis')
    plt.ylabel('Frequency')
    plt.legend()
    # plt.title('Radial Distribution of Tube Sprinkling')
    plt.tight_layout()
    plt.savefig(r'C:\Users\leehi\OneDrive\Documents\Imperial_tings\Fourth_Year\MSci Project\Thesis\Diagrams\Tubedistancer2+13+1distribution.png', dpi=300, transparent=True)
    plt.show()
