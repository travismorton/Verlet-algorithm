#!/usr/bin/python

import math
import numpy
import random

# initialize variables
dt = .01
boxWidth = 10.0
boxHeight = 10.0

#initialize arrays
r = numpy.array([[0, 0], [0, 1], [0, 2], [0, 3],
                [1, 0], [1, 1], [1, 2], [1, 3], [2, 0], [2, 1],
                [2, 2], [2, 3], [3, 0], [3, 1], [3, 2], [3, 3]])
r =r.astype(numpy.float64)
v = numpy.zeros((16, 2))
a = numpy.zeros((16, 2))

w = numpy.ndenumerate(r)
for x, y in w:
    r[x[0]][x[1]] += (random.random()-.5)/10.0

z = numpy.ndenumerate(v)
for x, y in z:
    v[x[0]][x[1]] = random.random()*3.0-1.5

print r
print v
#iterate through arrays with time step
