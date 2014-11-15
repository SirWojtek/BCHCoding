#!/usr/bin/python
from BinaryNumber import *
from copy import deepcopy

class CoderBCH:

	def __init__(self, generator, m):
		self._generator = generator
		self._m = m
		self._n = 2 ** m - 1;
		self._nk = self._generator.getPolyDegree()
		self._k = self._n - self._nk

	def encode(self, info):
		if not self._isEncodePossible(info):
			raise RuntimeError("Encoding imposible. Check length of message")
		return self._encode(info)

	def _isEncodePossible(self, info):
		return info.getPolyDegree() <= self._k

	def _encode(self, info):
		encMessage = deepcopy(info)
		encMessage.multiplyByX(self._nk)
		remainder = encMessage.divideRemainder(self._generator)
		encMessage.add(remainder)
		return encMessage

	def __repr__(self):
		return """Coder parameters:
generator:\t%s
m:\t\t%d
n:\t\t%d
k:\t\t%d
n - k:\t\t%d
""" % (self._generator, self._m, self._n, self._k, self._nk)

if __name__ == '__main__':
	gen = BinaryNumber(0b10011)
	info = BinaryNumber(0xA)
	coder = CoderBCH(gen, 4)
	print coder
	msg = coder.encode(info)
	print msg
