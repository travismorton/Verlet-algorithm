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

# iterate through arrays with time step
