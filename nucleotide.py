import matplotlib.pyplot as plt

def reverese_complement(seq):
    """Return the reverse complement of the input sequence."""
    seq = seq.upper()
    seq = seq.replace('A', 't')
    seq = seq.replace('T', 'a')
    seq = seq.replace('G', 'c')
    seq = seq.replace('C', 'g')
    return seq[::-1].upper()

def patternMatching(pattern,genome):
    positions = []
    for i in range(len(genome) - len(pattern) + 1):
        if genome[i:i+len(pattern)] == pattern:
            positions.append(i)
    return positions

def solveFromFile():
    with open('data/Vibrio_cholerae.txt') as f:
        genome = f.read()

    pattern = 'CTTGATCAT'
    positions = patternMatching(pattern, genome)
    for pos in positions:
        print(pos, end=' ')

def skew(genome):
    currentSkew = [0]
    for i in range(len(genome)):
        if genome[i] == 'C':
            currentSkew.append(currentSkew[i] - 1)
        elif genome[i] == 'G':
            currentSkew.append(currentSkew[i] + 1)
        else:
            currentSkew.append(currentSkew[i])
    return currentSkew

def minSkew(genome):
    skewList = skew(genome)
    minSkew = min(skewList)
    minSkewPos = []
    for i in range(len(skewList)):
        if skewList[i] == minSkew:
            minSkewPos.append(i)
    return minSkewPos

def hammingDistance(seq1, seq2):
    distance = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            distance += 1
    return distance

def approximatePatternMatching(pattern, genome, d):
    positions = []
    for i in range(len(genome) - len(pattern) + 1):
        if hammingDistance(pattern, genome[i:i+len(pattern)]) <= d:
            positions.append(i)
    return positions

def updatedFrequencyWords(genome, pattern, d):
    freqMap = {}
    n = len(genome)
    for i in range(n - len(pattern) + 1):
        if hammingDistance(pattern, genome[i:i+len(pattern)]) <= d:
            if pattern not in freqMap:
                freqMap[pattern] = 1
            else:
                freqMap[pattern] += 1
    return freqMap

def solveApproximatePatternMatching():
    with open('data/dataset_9_4.txt') as f:
        pattern = f.readline().strip()
        genome = f.readline().strip()
        d = int(f.readline().strip())
    positions = approximatePatternMatching(pattern, genome, d)
    for pos in positions:
        print(pos, end=' ')

def solveHammingDistance():
    with open('data/dataset_9_3.txt') as f:
        seq1 = f.readline().strip()
        seq2 = f.readline().strip()
    print(hammingDistance(seq1, seq2))

def plotSkewFromFile():
    with open('data/E_coli.txt') as f:
        genome = f.read()
    skewList = skew(genome)
    plt.plot(skewList)
    plt.show()

def solveMinSkew():
    with open('data/E_coli.txt') as f:
        genome = f.read()
    minSkews = minSkew(genome)
    for pos in minSkews:
        print(pos, end=' ')

# solveMinSkew()
# plotSkewFromFile()
# solveFromFile()
# solveHammingDistance()
# print(approximatePatternMatching('ATTCTGGA', 'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT', 3))
# solveApproximatePatternMatching()
print(updatedFrequencyWords('AACAAGCTGATAAACATTTAAAGAG', 'AAAAA', 1))