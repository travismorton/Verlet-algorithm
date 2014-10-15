#!/usr/bin/python

import math
import numpy
import random

# initialize variables
dt = .01
initialDx = .01
boxWidth = 10.0
boxHeight = 10.0
particleNumber = 16  # should be a square number
timeLength = 100

# initialize position, velocity,
# and acceleration arrays
r = numpy.zeros((particleNumber, 2))
v = numpy.zeros((particleNumber, 2))
a = numpy.zeros((particleNumber, 2))

# filling r array with equally spaced particles
w = numpy.ndenumerate(r)
count = 0
for x, y in w:
    xStep = boxWidth / float(math.sqrt(particleNumber))
    yStep = boxHeight / float(math.sqrt(particleNumber))
    if x[1] == 0:
        r[x[0]][x[1]] = count / int(math.sqrt(particleNumber)) * \
                        xStep + (random.random() - .5) * initialDx
        count += 1
    else:
        r[x[0]][x[1]] = x[0] % math.sqrt(particleNumber) * \
                        yStep + (random.random() - .5) * initialDx

# Fills velocity array with random velocities homogeneously
# distributed between -3/2 and 3/2
z = numpy.ndenumerate(v)
for x, y in z:
    v[x[0]][x[1]] = random.random() * 3.0 - 1.5

# finding distance between a list of points and one point
# and returns a list of distances
def distance(p0, p1, dim):
    return 0


# function for finding force/acceleration
def forcex(dx, dy):
    return -24 * dx * ((dx ** 2 + dy ** 2) ** 3 - 2) / (dx ** 2 + dy ** 2) ** 7


def forcey(dx, dy):
    return -24 * dy * ((dx ** 2 + dy ** 2) ** 3 - 2) / (dx ** 2 + dy ** 2) ** 7


# Define functions for iterating r, v, and a
def iterR(rcurrent, ):
    return 0


def iterV():
    return 0


def iterA():
    return 0

# iterate through arrays with time step
for t in range(int(timeLength / dt)):
    print t