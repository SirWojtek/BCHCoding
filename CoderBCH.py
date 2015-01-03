#!/usr/bin/python
from Polynomial import *
import numpy, copy

class CoderBCH:

    def __init__(self, fieldGenerator, generator, m, t):
        self._generator = generator
        self._m = m
        self._n = 2 ** m - 1;
        self._nk = self._generator.degree()
        self._k = self._n - self._nk
        self._t = t
        self._field = Field(m, fieldGenerator)

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
        # syndrome = info % self._generator
        # if syndrome.hammingWeight() == 0:
        #     print 'No transmission error'
        #     return info

        # assume that generator first minimal poly is m1
        t = self._getPartSyndroms(info , 1)
        print 'Part syndroms computed'

        fi = self._euclidian(t)
        print 'Euclidian algorithm ended'

        print self._getErrorPositions(fi)

    def _getErrorPositions(self, fi):
        result = []

        for i in range(self._n):
            if not Polynomial.getValueUsingAlphaMap(fi , i, self._field):
                result.append(self._field.getReversedPower(i))

        return result

    def _getPartSyndroms(self, info, m0):
        partSyndroms = {}
        polyAlpha = info.getAlphaMap()

        for i in range(2 * self._t):
            alphaPower = i + m0
            value = Polynomial.getValueUsingAlphaMap(polyAlpha, alphaPower, self._field)
            if value is None:
                print 'Syndrom S(alpha^%d) is 0' % alphaPower
            else:
                partSyndroms[i] = value

        return partSyndroms

    def _euclidian(self, T):
        old_r = Polynomial(1) * (2 * self._t)
        old_r = old_r.getAlphaMap()
        r = copy.deepcopy(T)
        old_t = {}
        t = {0 : 0}
        old_s = { 0 : 0}
        s = {}

        while Polynomial.getMapMaxKey(r) >= self._t:
            (result, remainder) = Polynomial.divideUsingAlphaMap(old_r, r, self._field)
            quotient = result

            old_r = r
            r = remainder

            temp = copy.deepcopy(old_s)
            old_s = s
            s = Polynomial.multiplyUsingAlphaMap(s, quotient, self._field)
            s = Polynomial.addUsingAlphaMap(s, temp, self._field)

            temp = copy.deepcopy(old_t)
            old_t = t
            t = Polynomial.multiplyUsingAlphaMap(t, quotient, self._field)
            t = Polynomial.addUsingAlphaMap(t, temp, self._field)

        # print "Bezout coefficients:"
        # print old_s, old_t
        # print "greatest common divisor:"
        # print old_r
        # print "quotients by the gcd:"
        # print t, s

        return t

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
    t = 12
    m = 8
    fieldGen = Polynomial(0435)
    gen = Polynomial(07500415510075602551574724514601)
    info = Polynomial(0b101011100110101110011010111001101011100110101110011010111001010)
    coder = CoderBCH(fieldGen, gen, m, t)
    print coder
    print 'INFO: ' + str(info)
    encodedMsg = coder.encode(info)
    print 'ENCODED: ' + str(encodedMsg)
    #Have some problem with generate proper noise.
    noisedMsg = addNoise(encodedMsg, t, gen.degree())
    #noisedMsg = encodedMsg
    print 'NOISED: ' + str(noisedMsg)
    decodedMsg = Polynomial()
    try:
        decodedMsg = coder.decode(noisedMsg)
    except RuntimeError, e:
        print e
    else:
        print 'DECODED: ' + str(decodedMsg)

    if decodedMsg == info:
        print 'INFO and DECODED messages match!'
    else:
        print 'No match at all'
    print '-----------------------------------------------'
    decodedMsgEuclid = coder.decodeEuclid(encodedMsg)
    print '-----------------------------------------------'
    print 'DECODED EUCLID: ' + str(decodedMsgEuclid)
