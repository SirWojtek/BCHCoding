#!/usr/bin/python
from numpy.polynomial import polynomial

class Polynomial:

	def __init__(self, number = 0):
		self._poly = self._getPolynomial(number)
		self._number = number

	def getPoly(self):
		return self._poly

	def getNumber(self):
		return self._number

	def divideRemainder(self, binaryNumber):
		resultPoly = polynomial.polydiv(self._poly, binaryNumber._poly)[1]
		resultPoly = self._normalizePoly(resultPoly)

		result = Polynomial()
		result._number = self._getNumberFromPoly(resultPoly)
		result._poly = resultPoly

		return result

	def getPolyDegree(self):
		return len(self._poly) - 1

	def multiplyByX(self, count):
		self._poly = [0] * count + self._poly
		self._number = self._getNumberFromPoly(self._poly)

	def add(self, binaryNumber):
		result = polynomial.polyadd(
			self._poly, binaryNumber._poly)
		result = self._normalizePoly(result)

		self._poly = result
		self._number = self._getNumberFromPoly(self._poly)


	def _getNumberFromPoly(self, poly):
		number = 0
		power = 0
		for digit in reversed(poly):
			number += 2 ** power * int(digit)
			power += 1
		return number

	def _normalizePoly(self, poly):
		norm = [abs(int(i)) % 2 for i in poly]
		return norm

	def _getPolynomial(self, number):
		reversedDigit = []
		for digit in reversed(self._translateToBinaryList(number)):
			reversedDigit.append(digit)
		return reversedDigit		

	def _translateToBinaryList(self, number):
		return [int(d) for d in str(bin(number))[2:]]

	def __repr__(self):
		return self._poly.__repr__()

if __name__ == '__main__':
	a = Polynomial(13)
	b = Polynomial(14)
	print a
	print b
	div = a.divideRemainder(b)
	print div

	a.add(b)
	print a
