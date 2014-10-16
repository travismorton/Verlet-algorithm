#!/usr/bin/python

import os
import math
import numpy
import random
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Initialize variables
dt = .001
initialDr = .01
boxWidth = 10.0
boxHeight = 10.0
particleNumber = 16  # should be a square number
timeLength = 2.5
vRange = 3.0
r2 = []  # average distance squared list, in order of time steps

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
                        xStep + (random.random() - .5) * initialDr
        count += 1
    else:
        r[x[0]][x[1]] = x[0] % math.sqrt(particleNumber) * \
                        yStep + (random.random() - .5) * initialDr

# Saving initial r array for use later in calculating
# total distance travelled, state, etc.
rInitial = numpy.copy(r)

# Fills velocity array with random velocities homogeneously
# distributed between -3/2 and 3/2
z = numpy.ndenumerate(v)
for x, y in z:
    v[x[0]][x[1]] = random.random() * 3.0 - 1.5


# functions for finding force/acceleration
def forcex(dx, dy):
    return -24 * dx * ((dx ** 2 + dy ** 2) ** 3 - 2) / (dx ** 2 + dy ** 2) ** 7


def forcey(dx, dy):
    return -24.0 * dy * ((dx ** 2 + dy ** 2) ** 3 - 2.0) / (dx ** 2 + dy ** 2) ** 7


# iterate through arrays with time step
for t in range(int(timeLength / dt)):
    # Adding pseudo-particles to each boundary of box to allow for
    # periodic forces
    rPseudo = numpy.copy(r)
    for particle in range(particleNumber):
        # left edge
        if vRange >= rPseudo[particle, 0] >= 0:
            rPseudo = numpy.append(rPseudo, [[r[particle, 0] + boxWidth, r[particle, 1]]], axis=0)
        # right edge
        if boxWidth >= rPseudo[particle, 0] >= boxWidth - vRange:
            rPseudo = numpy.append(rPseudo, [[r[particle, 0] - boxWidth, r[particle, 1]]], axis=0)
        # bottom edge
        if vRange >= rPseudo[particle, 1] >= 0:
            rPseudo = numpy.append(rPseudo, [[r[particle, 0], r[particle, 1] + boxHeight]], axis=0)
        # top edge
        if boxHeight >= rPseudo[particle, 1] >= boxHeight - vRange:
            rPseudo = numpy.append(rPseudo, [[r[particle, 0], r[particle, 1] - boxHeight]], axis=0)
        # bottom left corner
        if (vRange >= rPseudo[particle, 0] >= 0) and (vRange >= rPseudo[particle, 1] >= 0):
            rPseudo = numpy.append(rPseudo, [[r[particle, 0] + boxWidth, r[particle, 1] + boxHeight]], axis=0)
        # bottom right corner
        if (boxWidth >= rPseudo[particle, 0] >= boxWidth - vRange) and (vRange >= rPseudo[particle, 1] >= 0):
            rPseudo = numpy.append(rPseudo, [[r[particle, 0] - boxWidth, r[particle, 1] + boxHeight]], axis=0)
        # top left corner
        if (vRange >= rPseudo[particle, 0] >= 0) and (boxHeight >= rPseudo[particle, 1] >= boxHeight - vRange):
            rPseudo = numpy.append(rPseudo, [[r[particle, 0] + boxWidth, r[particle, 1] - boxHeight]], axis=0)
        # top right corner
        if (boxWidth >= rPseudo[particle, 0] >= boxWidth - vRange) and (boxHeight >= rPseudo[particle, 1] >=
                                                                                boxHeight - vRange):
            rPseudo = numpy.append(rPseudo, [[r[particle, 0] - boxWidth, r[particle, 1] - boxHeight]], axis=0)

    aPrev = numpy.copy(a)
    # Calculating force/acceleration matrix using potential function
    # and the pseudo particles created above for the boundary conditions
    for i in range(particleNumber):
        a[i][0] = 0
        a[i][1] = 0
        x = rPseudo[i][0]
        y = rPseudo[i][1]
        for j in range(rPseudo.shape[0]):
            # might need to switch particle i and j around...
            dx = x - rPseudo[j][0]
            dy = y - rPseudo[j][1]
            if math.sqrt(dx ** 2 + dy ** 2) < vRange and dx != 0 and dy != 0:
                a[i][0] += forcex(dx, dy)
                a[i][1] += forcey(dx, dy)

    vPrev = numpy.copy(v)
    # Calculating the velocities of each particle using the
    # current acceleration, previous acceleration, and previous velocity
    for i in range(particleNumber):
        v[i][0] = vPrev[i][0] + dt / 2 * (aPrev[i][0] + a[i][0])
        v[i][1] = vPrev[i][1] + dt / 2 * (aPrev[i][1] + a[i][1])

    rPrev = numpy.copy(r)
    # Calculating the r array using velocities and accelerations
    # found in the previous time step. Makes use of a periodic boundary
    # condition
    for i in range(particleNumber):
        r[i][0] = rPrev[i][0] + dt * vPrev[i][0] + dt ** 2 / 2 * aPrev[i][0]
        r[i][1] = rPrev[i][1] + dt * vPrev[i][1] + dt ** 2 / 2 * aPrev[i][1]
        if r[i][0] > boxWidth:
            r[i][0] -= boxWidth
        elif r[i][0] < 0:
            r[i][0] += boxWidth
        elif r[i][1] > boxHeight:
            r[i][1] -= boxHeight
        elif r[i][1] < 0:
            r[i][1] += boxHeight

    r2Placeholder = 0
    # Calculating distance (squared) each particle has travelled so far
    # and returning the average at this time step
    for i in range(particleNumber):
        dim = numpy.array([boxWidth, boxHeight])
        diff = numpy.abs(rInitial[i] - r[i])
        diff = numpy.where(diff > 0.5 * dim, dim - diff, diff)
        r2Placeholder += (diff[0] ** 2 + diff[1] ** 2) / particleNumber
    r2.append(r2Placeholder)

    p = 0
    # print a, v, and r to find t where spike is
    for i in range(particleNumber):
        if 1000 < abs(a[i, 0]):
            print i
            print aPrev
            print a
            p = 1
            break

    if p == 1:
        break
x = numpy.split(r, 2, axis=1)[0]
y = numpy.split(r, 2, axis=1)[1]

# Plotting (soon to come)
#def update_line(num, data, line):
#    line.set_data(data[..., :num])
#    return line,

print r2
#fig = plt.figure()
#l, = plt.scatter([], [])
#plt.xlim(0, 1)
#plt.ylim(0, 1)
#data = numpy.random.rand(2, 25)
#line_ani = animation.FuncAnimation(fig, update_line, 25,
#                                   fargs=(data, l), interval=50, blit=True)
#line_ani.save('lines.mp4')