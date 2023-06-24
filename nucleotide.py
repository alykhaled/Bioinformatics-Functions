import matplotlib.pyplot as plt
from frequentWords import frequentWords, betterFrequentWords, findClumps, frequencyTable,optimizedFindClumps

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
    count = 0
    for i in range(len(genome) - len(pattern) + 1):
        if hammingDistance(pattern, genome[i:i+len(pattern)]) <= d:
            count += 1
    return count

def immediateNeighbors(pattern):
    neighbors = []
    for i in range(len(pattern)):
        symbol = pattern[i]
        for nucleotide in ['A', 'C', 'G', 'T']:
            if nucleotide != symbol:
                neighbor = pattern[:i] + nucleotide + pattern[i+1:]
                neighbors.append(neighbor)
    return neighbors

def neighbors(pattern, d):
    if d == 0:
        return [pattern]
    if len(pattern) == 1:
        return ['A', 'C', 'G', 'T']
    neighborhood = []
    suffixNeighbors = neighbors(pattern[1:], d)
    for text in suffixNeighbors:
        if hammingDistance(pattern[1:], text) < d:
            for nucleotide in ['A', 'C', 'G', 'T']:
                neighborhood.append(nucleotide + text)
        else:
            neighborhood.append(pattern[0] + text)
    return neighborhood

def frequentWordsWithMismatches(genome, k, d):
    frequentPatterns = []
    neighborhoods = []
    for i in range(len(genome) - k + 1):
        neighborhoods.append(neighbors(genome[i:i+k], d))
    neighborhoods = [item for sublist in neighborhoods for item in sublist]
    neighborhoods = list(set(neighborhoods))
    counts = []
    for i in range(len(neighborhoods)):
        pattern = neighborhoods[i]
        counts.append(updatedFrequencyWords(genome, pattern, d))
    maxCount = max(counts)
    for i in range(len(neighborhoods)):
        if counts[i] == maxCount:
            frequentPatterns.append(neighborhoods[i])
    return frequentPatterns

def frequentWordsWithMismatchesAndReverseComplements(genome, k, d):
    frequentPatterns = []
    neighborhoods = []
    genome = genome + reverese_complement(genome)
    for i in range(len(genome) - k + 1):
        neighborhoods.append(neighbors(genome[i:i+k], d))
    neighborhoods = [item for sublist in neighborhoods for item in sublist]
    neighborhoods = list(set(neighborhoods))
    counts = []
    for i in range(len(neighborhoods)):
        pattern = neighborhoods[i]
        counts.append(updatedFrequencyWords(genome, pattern, d))
    maxCount = max(counts)
    for i in range(len(neighborhoods)):
        if counts[i] == maxCount:
            frequentPatterns.append(neighborhoods[i])
    return frequentPatterns

def solveFrequentWordsWithMismatchesAndReverseComplements():
    with open('data/dataset_9_10.txt') as f:
        genome = f.readline().strip()
        k, d = f.readline().strip().split()
    k = int(k)
    d = int(d)
    frequentPatterns = frequentWordsWithMismatchesAndReverseComplements(genome, k, d)
    for pattern in frequentPatterns:
        print(pattern, end=' ')

def solveFrequentWordsWithMismatches():
    with open('data/dataset_9_9.txt') as f:
        genome = f.readline().strip()
        k, d = f.readline().strip().split()
    k = int(k)
    d = int(d)
    frequentPatterns = frequentWordsWithMismatches(genome, k, d)
    for pattern in frequentPatterns:
        print(pattern, end=' ')

def solveNeighbors():
    with open('data/dataset_3014_4.txt') as f:
        pattern = f.readline().strip()
        d = int(f.readline().strip())
    neighborsList = neighbors(pattern, d)
    for neighbor in neighborsList:
        print(neighbor, end=' ')

def solveUpdatedFrequencyWords():
    with open('data/dataset_9_6.txt') as f:
        pattern = f.readline().strip()
        genome = f.readline().strip()
        d = int(f.readline().strip())
    print(updatedFrequencyWords(genome, pattern, d))
    

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

def readFasta(filename):
    with open(filename) as f:
        genome = ''
        for line in f:
            if line[0] != '>':
                genome += line.strip()
    return genome

def solveSalmonella():
    genome = readFasta('data/Salmonella_enterica.txt')
    minSkews = minSkew(genome)
    minGenomes = genome[minSkews[0]:minSkews[0]+500]
    print(len(minGenomes))
    boxes = frequentWordsWithMismatchesAndReverseComplements(minGenomes, 9, 1)
    for box in boxes:
        print(box, end=' ')

# solveMinSkew()
# plotSkewFromFile()
# solveFromFile()
# solveHammingDistance()
# print(approximatePatternMatching('ATTCTGGA', 'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT', 3))
# solveApproximatePatternMatching()
# print(updatedFrequencyWords('TTTAGAGCCTTCAGAGG', 'GAGG', 2))
# solveUpdatedFrequencyWords()
# print(immediateNeighbors('ACG'))
# print(neighbors('ACG', 1))
# solveNeighbors()
# print(frequentWordsWithMismatches('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4, 1))
# solveFrequentWordsWithMismatches()
# print(frequentWordsWithMismatchesAndReverseComplements('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4, 1))
# solveFrequentWordsWithMismatchesAndReverseComplements()
solveSalmonella()