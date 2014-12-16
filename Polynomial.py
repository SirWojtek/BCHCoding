#!/usr/bin/python
from numpy import roll
from numpy.polynomial import polynomial

class Field:
    def __init__(self, degree = 0, generator = None):
        if generator is not None:
            self._degree = generator.degree()
            self._generator = generator
        else:
            self._degree = degree
            self._generator = self._createGenerator(degree)
        self._alphaPowerMap = self._getPowerMap()

    def getAlpha(self, degree):
        if degree < 0:
            degree = pow(2, self._degree) - 1 - degree

        poly = Polynomial(0b1)*degree
        return poly % self._generator

    def getAlphaPower(self, poly):
        return self._alphaPowerMap.get(poly, None)

    def getGenerator(self):
        return self._generator

    # W razie potrzeby wiekszych: http://theory.cs.uvic.ca/gen/poly.html
    def _createGenerator(self, fieldDegree):
        minimalPolys = [0b11, 0b111, 0b1011, 0b10011, 0b100101, 0b1000011,
                        0b10000011, 0b101101001, 0b1000000011, 0b10000001001]
        try:
            return Polynomial(minimalPolys[fieldDegree-1])
        except IndexError:
            raise ValueError('There is no minimal polynomial in table for field degree = ' + str(fieldDegree))

    def _getPowerMap(self):
        result = {}
        for i in range(pow(2, self._degree) - 1):
            result[self.getAlpha(i)] = i
        return result

class Polynomial:

    def __init__(self, number = 0, poly = None):
        if poly is not None:
            self._poly = poly
        else:
            self._poly = self._getPolynomial(number)

    def degree(self):
        try:
            return (len(self._poly) - 1) - self._poly[::-1].index(1)
        except ValueError:
            return 0

    def hammingWeight(self):
        weight = 0
        for bit in self._poly:
            if bit > 0:
                weight += 1
        return weight

    def copy(self):
        return Polynomial(poly=list(self._poly))

    # function returns poly(alpha^alphaPower)
    def polyValAlpha(self, alphaPower, m):
        maxAlphaPower = pow(2, m) - 2
        self._trimPoly(self._poly)
        poly = self._substitiuteWithAlphaPower(alphaPower, maxAlphaPower)
        return self._sumAlpha(poly, m)

    def _substitiuteWithAlphaPower(self, alphaPower, maxAlphaPower):
        maxDegree = self._getMaxDegree(alphaPower, maxAlphaPower)
        result = [0 for i in range(maxDegree + 1)]

        for degree in range(len(self._poly)):
            currentIndex = degree
            newIndex = (alphaPower * degree) % (maxAlphaPower + 1)

            result[newIndex] += self._poly[currentIndex]
        result = self._normalizePoly(result, 0)
        self._trimPoly(result)
        return result

    def _getMaxDegree(self, alphaPower, maxAlphaPower):
        if alphaPower * self.degree() > maxAlphaPower:
            return maxAlphaPower
        else:
            return alphaPower * self.degree()

    def _sumAlpha(self, poly, m):
        field = Field(m)
        result = Polynomial(0)

        for i in range(len(poly)):
            if poly[i] == 1:
                result += field.getAlpha(i)

        return field.getAlphaPower(result)

    def divideByAlphaMap(self, division, m):
        field = Field(m)
        divisor = self._getAlphaMapDivisor()
        return self._divideUsingAlphaMap(field, divisor, division)

    def _getAlphaMapDivisor(self):
        result = {}
        for i in range(len(self._poly)):
            if self._poly[i] == 1:
                result[i] = 0
            else:
                result[i] = None
        return result

    def _divideUsingAlphaMap(self, field, divisor, division):
        # TODO: implement this awesome algorithm
        return None

    def _normalizePoly(self, poly, expectedSize):
        norm = [abs(int(i)) % 2 for i in poly]
        sizeDiff = expectedSize - len(norm)
        if sizeDiff > 0:
            norm += [0]*sizeDiff
        return norm

    def _trimPoly(self, poly):
        while len(poly) and poly[-1] == 0:
            poly.pop()

    def _getPolynomial(self, number):
        reversedDigit = []
        for digit in reversed(self._translateToBinaryList(number)):
            reversedDigit.append(digit)
        return reversedDigit

    def _getNumberFromPoly(self):
        number = 0
        power = 0
        for digit in self._poly:
            number += 2 ** power * int(digit)
            power += 1
        return number

    def _translateToBinaryList(self, number):
        return [int(d) for d in str(bin(number))[2:]]

    def __len__(self):
        return len(self._poly)

    def __getitem__(self, key):
        return self._poly[key]

    def __setitem__(self, key, value):
        self._poly[key] = value

    def __str__(self):
        return ' '.join(list(str(bin(self._getNumberFromPoly()))[2:]))

    def __repr__(self):
        return self._poly.__repr__()

    def __int__(self):
        return self._getNumberFromPoly()

    def __hex__(self):
        return hex(self._getNumberFromPoly())

    def __oct__(self):
        return oct(self._getNumberFromPoly())

    def __index__(self):
        return self._getNumberFromPoly()

    def __mod__(self, other):
        result = polynomial.polydiv(self._poly, other._poly)[1]
        result = self._normalizePoly(result, len(self._poly))
        return Polynomial(poly=result)

    def __mul__(self, other):
        if isinstance(other, Polynomial):
            result = polynomial.polymul(self._poly, other._poly)
        else:
            result = [0] * other + self._poly
        return Polynomial(poly=result)

    def __imul__(self, other):
        return self.__mul__(other)

    def __div__(self, other):
        if isinstance(other, Polynomial):
            result = polynomial.polydiv(self._poly, other._poly)[0]
        else:
            result = self._poly[other:]
        result = self._normalizePoly(result, len(self._poly))
        return Polynomial(poly=result)

    def __idiv__(self, other):
        return self.__div__(other)

    def __add__(self, other):
        result = polynomial.polyadd(self._poly, other._poly)
        result = self._normalizePoly(result, len(self._poly))
        return Polynomial(poly=result)

    def __iadd__(self, other):
        return self.__add__(other)

    def __lshift__(self, other):
        result = list(roll(self._poly, other))
        return Polynomial(poly=result)

    # and this function too
    def __rshift__(self, other):
        result = list(roll(self._poly, -other))
        return Polynomial(poly=result)

    def __eq__(self, other):
        return self._getNumberFromPoly() == other._getNumberFromPoly()

    def __hash__(self):
        return self._getNumberFromPoly()

if __name__ == '__main__':
    a = Polynomial(0b1000)
    b = Polynomial(0b0011)
    c = Polynomial(0)
    print a
    print b
    print a*8
    print a+b
    print a/2
    print a%b
    print a*b

    # print '----------------------'


    # f = Field(8)
    # f.getAlphaPower(c)
