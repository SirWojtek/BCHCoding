#!/usr/bin/python
from Polynomial import *
import numpy, copy

class CoderBCH:

    def __init__(self, generator, m, t):
        self._generator = generator
        self._m = m
        self._n = 2 ** m - 1;
        self._nk = self._generator.degree()
        self._k = self._n - self._nk
        self._t = t
        self._field = Field(m)

    def encode(self, info):
        if not self._isEncodePossible(info):
            raise RuntimeError("Encoding imposible. Check length of message")
        return self._encode(info)

    def decode(self, info):
        return self._decode(info) / self._nk

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
        print 'Part syndroms computed'

        B = self._euclidian(t)
        print 'Euclidian algorithm ended'

        delta = {0 : Polynomial.getReversedPower(self._m, B[0]) }
        (fi, rem) = Polynomial.divideUsingAlphaMap(self._m, B, delta)
        print 'Fi function computed'

        print self._getErrorPositions(fi)

    def _getErrorPositions(self, fi):
        result = []

        for i in range(self._n):
            if not Polynomial.getValueUsingAlphaMap(self._m, fi , i):
                result.append(i)

        return result


    def _getPartSyndroms(self, syndrome, m0):
        partSyndroms = {}
        for i in range(m0, 2 * self._t - 1):
            polyAlpha = syndrome.polyValAlpha(i, self._m)
            if polyAlpha:
                partSyndroms[i] = polyAlpha
            else:
                print 'Syndrom S%d is 0' % i
        return partSyndroms

    def _euclidian(self, T):
        A = [[ 1, 0 ],
            [ 0, 1 ]]
        s = Polynomial(1) * (2 * self._t)
        s = s.getAlphaMap()
        t = T

        while Polynomial.getMapMaxKey(t) >= self._t:
            AOld = copy.deepcopy(A)

            (result, remainder) = Polynomial.divideUsingAlphaMap(self._m, s, t)
            Q = result
            s = t
            t = remainder

            A[0][0] = AOld[1][0]
            A[0][1] = AOld[1][1]
            minusQ = self._getMinusQ(Q)
            # A[0][0] + A[1][0] * Q
            A[1][0] = self._getA(AOld[0][0], AOld[1][0], minusQ)
            # A[0][1] + A[1][1] * Q
            A[1][1] = self._getA(AOld[0][1], AOld[1][1], minusQ)

        return A[1][1]

    def _getMinusQ(self, Q):
        polyAlphaZero = {}
        return Polynomial.subUsingAlphaMap(self._m, polyAlphaZero, Q)

    def _getA(self, a0, a1, Q):
        if a0 == 0 and a1 == 1:
            return Q
        if a0 == 1 and a1 == 0:
            return 1
        if a0 == 1:  # a1 == alphaMap
            resultQ = Polynomial.multiplyUsingAlphaMap(self._m, Q, a1)
            return Polynomial.addUsingAlphaMap(self._m, resultQ, {0 : 0}) # <- 1
        else: # a0 == alphaMap and a1 == alphaMap
            resultQ = Polynomial.multiplyUsingAlphaMap(self._m, Q, a1)
            return Polynomial.addUsingAlphaMap(self._m, resultQ, a0)



    def __repr__(self):
        return """Coder parameters:
generator:\t%s
m:\t\t%d
n:\t\t%d
k:\t\t%d
n - k:\t\t%d
t:\t\t%d
""" % (int(self._generator), self._m, self._n, self._k, self._nk, self._t)

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
    print '-----------------------------------------------'
    decodedMsgEuclid = coder.decodeEuclid(noisedMsg)
    print '-----------------------------------------------'
    print 'DECODED EUCLID: ' + str(decodedMsgEuclid)
