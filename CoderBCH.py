#!/usr/bin/python
from Polynomial import *
import numpy

class CoderBCH:

    def __init__(self, generator, m, t):
        self._generator = generator
        self._m = m
        self._n = 2 ** m - 1;
        self._nk = self._generator.degree()
        self._k = self._n - self._nk
        self._t = t

    def encode(self, info):
        if not self._isEncodePossible(info):
            raise RuntimeError("Encoding imposible. Check length of message")
        return self._encode(info)

    def decode(self, info):
        return self._decode(info) >> self._nk

    def _isEncodePossible(self, info):
        return info.degree() <= self._k

    def _encode(self, info):
        encodedMsg = info*self._nk
        return encodedMsg + (encodedMsg % self._generator)

    def _decode(self, info):
        decodedMsg = info.copy()
        for i in range(self._k):
            syndrome = decodedMsg % self._generator
            if syndrome.hammingWeight() <= self._t:
                print 'DECODE ITERATIONS: ' + str(i)
                decodedMsg = decodedMsg + syndrome
                decodedMsg = decodedMsg << i
                return decodedMsg
            decodedMsg = decodedMsg >> 1
        raise RuntimeError('Unable to correct errors for input message: ' + str(hex(info)))

    def decodeEuclid(self, info):
        syndrome = info % self._generator
        if syndrome.hammingWeight() == 0:
            print 'No transmission error'
            return info

        # assume that generator first minimal poly is m1
        t = self._getPartSyndroms(syndrome , 1)
        A = [[ 1, 0 ],
            [ 0, 1 ]]

        self._euclidian(t, A)

        # delta = A[1][1](0)  #  Bxx(0) from Euclidian algoritm (recurse 0)
        # A = partSyndroms * delta ^ -1
        # fi = A[1][1] * delta ^ -1


    def _getPartSyndroms(self, syndrome, m0):
        partSyndroms = []
        for i in range(m0, 2 * t - 1):
            polyAlpha = syndrome.polyValAlpha(i, pow(2, m) - 2)
            if polyAlpha.degree() != 0:
                partSyndroms.append(syndrome.polyValAlpha(i, pow(2, m) - 2))
            else:
                print 'Syndrom S%d is 0' % i
        return partSyndroms

    def _euclidian(self, t, A):
        s = Polynomial(1) << (2 * self._t)

        # while len(t) - 1 >= self._t:
        #     sOld = s
        #     tOld = t
        #     AOld = A

        #     Q = sOld / tOld
        #     s = tOld
        #     t = sOld + Q * tOld

        #     A[0][0] = AOld[1][0]
        #     A[0][1] = AOld[1][1]
        #     A[1][0] = AOld[0][0] + AOld[1][0] * Q
        #     A[1][1] = AOld[0][1] +  AOld[1][1] * Q

    def __repr__(self):
        return """Coder parameters:
generator:\t%s
m:\t\t%d
n:\t\t%d
k:\t\t%d
n - k:\t\t%d
t:\t\t%d
""" % (int(self._generator), self._m, self._n, self._k, self._nk, self._t)


def addNoise(message, maxErrors, generatorDegree=None):
    poly = message.copy()
    if generatorDegree is not None:
        positions = numpy.random.randint(generatorDegree, size=maxErrors)
    else:
        positions = numpy.random.randint(len(poly), size=maxErrors)
    positions = sorted(list(set(positions)))
    print 'NOISE_POSITIONS' + str(positions)
    for position in positions:
        poly[position] = int(not(poly[position]))
    return poly


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
    noisedMsg = addNoise(encodedMsg, t, gen.degree())
    #noisedMsg = encodedMsg
    print 'NOISED: ' + str(noisedMsg)
    decodedMsg = coder.decode(noisedMsg)
    print 'DECODED: ' + str(decodedMsg)
    if decodedMsg == info:
        print 'INFO and DECODED messages match!'
    else:
        print 'No match at all'
    decodedMsgEuclid = coder.decodeEuclid(noisedMsg)
    print 'DECODED EUCLID: ' + str(decodedMsgEuclid)
