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

solveFromFile()