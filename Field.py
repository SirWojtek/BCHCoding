#!/usr/bin/python
from Polynomial import *

class Field:

    def __init__(self, degree = 0, generator = None):
        if generator is not None:
            self._degree = generator.degree()
            self._generator = generator
        else:
            self._degree = degree
            self._generator = self._createGenerator(degree)

    def getAlpha(self, degree):
        poly = Polynomial(0b1)*degree
        return poly % self._generator

    def getGenerator(self):
        return self._generator

    # W razie potrzeby wiekszych: http://theory.cs.uvic.ca/gen/poly.html
    def _createGenerator(self, fieldDegree):
        minimalPolys = [0b11, 0b111, 0b1011, 0b10011, 0b100101, 0b1000011,
                        0b10000011, 0b100011011, 0b1000000011, 0b10000001001]
        try:
            return Polynomial(minimalPolys[fieldDegree-1])
        except IndexError:
            raise ValueError('There is no minimal polynomial in table for field degree = ' + str(fieldDegree))


if __name__ == '__main__':
    f = Field(4)
    print f.getAlpha(16)
    print f.getGenerator()
