from CoderBCH import CoderBCH
from Polynomial import Polynomial
import numpy

codes = [
    {
        't' : 12,
        'm' : 8,
        'fieldGen' : Polynomial(0435),
        'codeGen' : Polynomial(07500415510075602551574724514601),
        'info' : Polynomial(0b1010111001101011100110101110010101011101011100111011011011010100101011111101101010111001101011100101001010111100110101010101111001101011100110101110011010111001010)
    },
    {
        't' : 55,
        'm' : 8,
        'fieldGen' : Polynomial(0435),
        'codeGen' : Polynomial(01256215257060332656001773153607612103227341405653074542521153121614466513473725),
        'info' : Polynomial(0b101011110111111011010)
    },
    {
        't' : 29,
        'm' : 8,
        'fieldGen' : Polynomial(0435),
        'codeGen' : Polynomial(024024710520644321515554172112331163205444250362557643221706035),
        'info' : Polynomial(0b10101110010100101011110011010101010111100110101110011010111001101011100)
    }
]

def randErrorPositions(range, numberOfErrors):
    if range < numberOfErrors:
        raise RuntimeError("Range value cant be bigger then number of errors!")
    positions = list()
    diff = numberOfErrors
    while diff > 0:
        positions += list(numpy.random.randint(range, size=diff))
        positions = list(set(positions))
        diff = numberOfErrors - len(positions)
    return sorted(positions)

def moveRagne(positions, messageLen, errorsDelta):
    temp = list(positions)
    if len(positions) == 1:
        return temp
    diff = temp[-1] - temp[0]
    while diff <= errorsDelta:
        if temp[-1] < (messageLen-1):
            temp[-1] += 1
        if temp[0] > 0:
            temp[0] -= 1
        diff = temp[-1] - temp[0]
    return temp

def addNoise(message, numberOfErrors, errorsDelta, higherDelta=None):
    poly = message.copy()
    positions = randErrorPositions(errorsDelta, numberOfErrors)
    shift = numpy.random.randint(len(message) - errorsDelta)
    for i, pos in enumerate(positions):
        positions[i]  = pos + shift
    if higherDelta is not None:
        positions = moveRagne(positions, len(message), errorsDelta)
    for position in positions:
        poly[position] = int(not(poly[position]))
    return poly

def toCsv(path, result):
    file = open(path, 'w+')
    for t in result:
        file.write(str(t))
        for value in result[t]:
            file.write(','+ str(value))
        file.write('\n')
    file.close()

def testWithChangingT(code, numberOfTests=1):
    results = dict()
    coder = CoderBCH(code['fieldGen'], code['codeGen'], code['m'], code['t'])
    encodedMsg = coder.encode(code['info'])
    for t in range(1,code['t'] + 2):
        results[t] = list()
        for i in range(numberOfTests):
            #nosiedMsg = addNoise(encodedMsg, t, coder._nk, higherDelta=True)
            nosiedMsg = addNoise(encodedMsg, t, coder._nk)
            try:
                decodedMsg = coder.decode(nosiedMsg)
            except RuntimeError:
                results[t].append(0)
            else:
                results[t].append(int(decodedMsg == code['info']))
    return results

if __name__ == '__main__':
    results = testWithChangingT(codes[0], 100)
    toCsv('result.csv', results)
    results = testWithChangingT(codes[1], 100)
    toCsv('result2.csv', results)
    results = testWithChangingT(codes[2], 100)
    toCsv('result3.csv', results)