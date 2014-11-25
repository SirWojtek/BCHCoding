#!/usr/bin/python
from Polynomial import *
from copy import deepcopy
import numpy

class CoderBCH:

    def __init__(self, generator, m, t):
        self._generator = generator
        self._m = m
        self._n = 2 ** m - 1;
        self._nk = self._generator.getPolyDegree()
        self._k = self._n - self._nk
        self._t = t

    def encode(self, info):
        if not self._isEncodePossible(info):
            raise RuntimeError("Encoding imposible. Check length of message")
        return self._encode(info)

    def decode(self, info):
        msg = self._decode(info)
        msg.divideByX(self._nk)
        return msg

    def _isEncodePossible(self, info):
        return info.getPolyDegree() <= self._k

    def _encode(self, info):
        encMessage = deepcopy(info)
        encMessage.multiplyByX(self._nk)
        remainder = encMessage.divideRemainder(self._generator)
        encMessage.add(remainder)
        return encMessage

    def _decode(self, info):
        decMessage = deepcopy(info)
        for i in range(self._k):
            syndrome = decMessage.divideRemainder(self._generator)
            hammingWeight = self._getHammingWeight(syndrome)
            if hammingWeight <= self._t:
                print 'DECODE ITERATIONS: ' + str(i)
                decMessage.add(syndrome)
                decMessage.shiftLeft(i)
                return decMessage
            decMessage.shiftRight()
        raise TypeError('Unable to correct errors for input message: ' + str(bin(info.getNumber())))


    def _getHammingWeight(self, polynomial):
        weight = 0
        for bit in polynomial.getPoly():
            if bit > 0:
                weight += 1
        return weight

    def __repr__(self):
        return """Coder parameters:
generator:\t%s
m:\t\t%d
n:\t\t%d
k:\t\t%d
n - k:\t\t%d
t:\t\t%d
""" % (self._generator, self._m, self._n, self._k, self._nk, self._t)

def addNoise(message, maxErrors):
    poly = message.getPoly()
    positions = numpy.random.randint(len(poly), size=maxErrors)
    positions = list(set(positions))
    for position in positions:
        poly[position] = int(not(poly[position]))
    return Polynomial(poly=poly)


if __name__ == '__main__':
    t = 30
    m = 8
    gen = Polynomial(010754475055163544325315217357707003666111726455267613656702543301)
    info = Polynomial(0b101011100110101110011010111001101011100110101110011010111001010)
    coder = CoderBCH(gen, m, t)
    print coder
    print 'INFO: ' + str(info)
    encodedMsg = coder.encode(info)
    print 'ENCODED: ' + str(encodedMsg)
    #Have some problem with generate proper noise. 
    noisedMsg = addNoise(encodedMsg, 4)
    #noisedMsg = encodedMsg
    print 'NOISED: ' + str(noisedMsg)
    decodedMsg = coder.decode(noisedMsg)
    print 'DECODED: ' + str(decodedMsg)
    if decodedMsg.getNumber() == info.getNumber():
        print 'INFO and DECODED messages match!'
    else:
        print 'No match at all'
