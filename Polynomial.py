#!/usr/bin/python
from numpy import roll
from numpy.polynomial import polynomial

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
    # computed alpha power is reduced mod maxAlphaPower
    def polyValAlpha(self, alphaPower, maxAlphaPower):
        maxDegree = self._getMaxDegree(alphaPower, maxAlphaPower)
        result = [0 for i in range(maxDegree + 1)]
        for degree in range(len(self._poly)):
            currentIndex = degree
            newIndex = (alphaPower * degree) % (maxAlphaPower + 1)
            result[newIndex] += self._poly[currentIndex]
        result = self._normalizePoly(result, 0)
        return Polynomial(poly = result)

    def _getMaxDegree(self, alphaPower, maxAlphaPower):
        if alphaPower * self.degree() > maxAlphaPower:
            return maxAlphaPower
        else:
            return alphaPower * self.degree()

    def _normalizePoly(self, poly, expectedSize):
        norm = [abs(int(i)) % 2 for i in poly]
        sizeDiff = expectedSize - len(norm)
        if sizeDiff > 0:
            norm += [0]*sizeDiff
        return norm

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
        #TODO: write mul for polynomials
        if isinstance(other, Polynomial):
            result = [0] * other + self._poly
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

if __name__ == '__main__':
    a = Polynomial(0b1000)
    b = Polynomial(0b0011)
    print a+a
    '''
    print a
    print b
    print a*8
    print a+b
    print a/2
    print a%b
    print c.polyValAlpha(5, 14)
    '''
