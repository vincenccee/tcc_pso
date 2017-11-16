import matplotlib                  ## Bibliotecas de plot
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import time
import csv

import random
import math
import numpy as np

if __name__ == "__main__":

    CHANGES = 1

    # get the contour map
    reader = csv.reader(open("contour_data_1.csv", "r"), delimiter=",")
    x = list(reader)
    result = np.array(x).astype("float")

    # get the population particles
    reader = csv.reader(open("population_2.csv", "r"), delimiter=",")
    x = list(reader)
    population = np.array(x).astype("float")

    h, w = result.shape
    lb = 0
    ub = 100

    xlist = np.linspace(lb, ub, num=h)
    ylist = np.linspace(lb, ub, num=h)
    # levels = np.linspace(40, 55, num=50)

    colors = ['#ff0000', '#ff6600', '#ffff00', '#13f35d', '#59b578',
    '#0e7d33', '#008cff', '#002aff', '#1b1059', '#444444']
    cmap = LinearSegmentedColormap.from_list('mycmap', list(reversed(colors)))

    plt.figure()
    cp = plt.contourf(xlist, ylist, result, 50, cmap=cmap)
    plt.colorbar(cp)
    plt.plot(population[:,1], population[:,0], ".", color="black")
    plt.savefig("contuor_initial.png")

    # for i in range(0, CHANGES):
    #     filename1 = "population_"+str(i)+".csv"
    #     reader = csv.reader(open(filename1, "r"), delimiter=",")
    #     x = list(reader)
    #     result = np.array(x).astype("float")
    #
    #     filename2 = "contour_data_"+str(i+1)+".csv"
    #     reader = csv.reader(open(filename2, "r"), delimiter=",")
    #     x = list(reader)
    #     population = np.array(x).astype("float")
    #
    #     plt.figure()
    #     cp = plt.contourf(xlist, ylist, result, 50, cmap=cmap)
    #     plt.colorbar(cp)
    #     plt.plot(population[:,1], population[:,0], ".", color="black")
    #     filename1 = "contour_"+str(i)+"_change.csv"
    #     plt.savefig("contuor_initial.png")
    #
    # # get the contour map
    # reader = csv.reader(open("contour_data_"+str(CHANGES)+".csv", "r"), delimiter=",")
    # x = list(reader)
    # result = np.array(x).astype("float")
    #
    # # get the population particles
    # reader = csv.reader(open("population_"+str(CHANGES+1)+".csv", "r"), delimiter=",")
    # x = list(reader)
    # population = np.array(x).astype("float")
