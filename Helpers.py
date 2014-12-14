#!/usr/bin/python
from Polynomial import *
import numpy

#@param fieldDegree - stopien ciala GF(2^fieldDegree)
def getFieldGenerator(fieldDegree):
   minimalPolys = [0b11, 0b111, 0b1011, 0b10011, 0b100101, 0b1000011,
                   0b10001001, 0b100011101]
   try:
        return Polynomial(minimalPolys[fieldDegree-1])
   except IndexError:
       raise ValueError('There is no minimal polynomial in table for field degree = ' + str(fieldDegree))



def addNoise(message, maxErrors, generatorDegree):
    poly = message.copy()
    positions = numpy.random.randint(generatorDegree, size=maxErrors)
    shift = numpy.random.randint(len(message) - generatorDegree)
    positions = sorted(list(set(positions)))
    print len(positions)
    for i, pos in enumerate(positions):
        positions[i]  = pos + shift
    print 'NOISE_POSITIONS' + str(positions)
    for position in positions:
        poly[position] = int(not(poly[position]))
    return poly

if __name__ == '__main__':
    print getFieldGenerator(9)