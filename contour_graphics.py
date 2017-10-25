import matplotlib                  ## Bibliotecas de plot
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import time
import csv


import random
import math
import numpy as np

if __name__ == "__main__":

    reader = csv.reader(open("contour_data_0.csv", "r"), delimiter=",")
    x = list(reader)
    result = np.array(x).astype("float")

    h, w = result.shape
    lb = 0
    ub = 100

    xlist = np.linspace(lb, ub, num=h)
    ylist = np.linspace(lb, ub, num=h)

    colors = ['#ff0000', '#ff6600', '#ffff00', '#43ff00', '#00ff1d',
              '#00ffe1', '#008cff', '#002aff', '#1b1059', '#444444']
    cmap = LinearSegmentedColormap.from_list('mycmap', list(reversed(colors)))

    plt.figure()
    cp = plt.contourf(xlist, ylist, result, 50, cmap=cmap)
    plt.colorbar(cp)
    # plt.plot([a for a, b in p], [b for a, b in p], ".", color="black")
    plt.show()
    plt.savefig("figure.png")
