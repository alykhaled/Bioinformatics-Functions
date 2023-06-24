from count import patternCount

def frequentWords(text, k):
    """
    Finds the most frequent k-mers in a string
    :param text: the text to search
    :param k: the length of the k-mer
    :return: a list of the most frequent k-mers in the text
    """

    frequentPatterns = []
    count = []
    for i in range(len(text) - k + 1):
        pattern = text[i:i+k]
        count.append(patternCount(text, pattern))
    maxCount = max(count)
    for i in range(len(text) - k + 1):
        if count[i] == maxCount:
            frequentPatterns.append(text[i:i+k])
    return list(set(frequentPatterns))

def frequencyTable(text,k):
    """
    Creates a frequency table for a given text
    :param text: the text to search
    :param k: the length of the k-mer
    :return: a list of the number of times each k-mer appears in the text
    """
    freqMap = {}
    n = len(text)
    for i in range(n-k+1):
        pattern = text[i:i+k]
        if pattern not in freqMap:
            freqMap[pattern] = 1
        else:
            freqMap[pattern] += 1
    return freqMap

def betterFrequentWords(text, k):
    """
    Finds the most frequent k-mers in a string
    :param text: the text to search
    :param k: the length of the k-mer
    :return: a list of the most frequent k-mers in the text
    """
    frequentPatterns = []
    freqMap = frequencyTable(text, k)
    maxFreq = max(freqMap.values())
    for pattern in freqMap:
        if freqMap[pattern] == maxFreq:
            frequentPatterns.append(pattern)
    return frequentPatterns

def findClumps(genome, k, L, t):
    """
    Finds patterns forming clumps in a string
    :param genome: the text to search
    :param k: the length of the k-mer
    :param L: the length of the window
    :param t: the minimum number of times a k-mer must appear in the window
    :return: a list of the most frequent k-mers in the text
    """
    frequentPatterns = []
    n = len(genome)
    for i in range(n-L+1):
        window = genome[i:i+L]
        freqMap = frequencyTable(window, k)
        for pattern in freqMap:
            if freqMap[pattern] >= t:
                frequentPatterns.append(pattern)
    return list(set(frequentPatterns))

def optimizedFindClumps(genome, k, L, t):
    """
    Finds patterns forming clumps in a string
    :param genome: the text to search
    :param k: the length of the k-mer
    :param L: the length of the window
    :param t: the minimum number of times a k-mer must appear in the window
    :return: a list of the most frequent k-mers in the text
    """
    allPossibleKmers = []
    kmersIndex = {}
    for i in range(len(genome) - k + 1):
        kmer = genome[i:i+k]
        if kmer not in kmersIndex:
            kmersIndex[kmer] = [i]
        else:
            kmersIndex[kmer].append(i)
    for kmer in kmersIndex:
        if len(kmersIndex[kmer]) >= t:
            allPossibleKmers.append(kmer)
    frequentPatterns = []
    for kmer in allPossibleKmers:
        for i in range(len(kmersIndex[kmer]) - t + 1):
            if kmersIndex[kmer][i+t-1] - kmersIndex[kmer][i] <= L - k:
                frequentPatterns.append(kmer)
                break
    return list(set(frequentPatterns))

def solveFromFile():
    with open('data/E_coli.txt') as f:
        genome = f.read()
        k = 9
        L = 500
        t = 3
        # print(len(genome))
        clumps = optimizedFindClumps(genome, k, L, t)
        print(' '.join(clumps))
        print(len(clumps))


# print(frequentWords("ACGTTGCATGTCGCATGATGCATGAGAGCT", 3))
# print(frequencyTable("ACGTTGCATGTCGCATGATGCATGAGAGCT", 3))
# print(betterFrequentWords("CCTATACACCTATACATTGCGTAACCTATACACCTATACATCCGCCGTATCGTCTATCGTCTTGCGTAACCTATACATCCGCCGTCCGCCGGTAAACGTCCGCCGTCCGCCGTCCGCCGTATCGTCCCTATACATATCGTCTTGCGTAATCCGCCGCCTATACATATCGTCGTAAACGTCCGCCGGTAAACGTTGCGTAACCTATACAGTAAACGTATCGTCTCCGCCGTATCGTCTTGCGTAATATCGTCCCTATACAGTAAACGTATCGTCTCCGCCGCCTATACATATCGTCCCTATACAGTAAACGGTAAACGTATCGTCTCCGCCGTTGCGTAACCTATACATATCGTCTTGCGTAATATCGTCCCTATACATTGCGTAATCCGCCGTATCGTCTTGCGTAACCTATACATATCGTCCCTATACAGTAAACGGTAAACGTTGCGTAACCTATACATATCGTCGTAAACGTCCGCCGTCCGCCGCCTATACACCTATACATCCGCCGCCTATACAGTAAACGGTAAACGCCTATACATATCGTCTCCGCCGCCTATACATATCGTCTCCGCCGCCTATACATATCGTCTTGCGTAATATCGTCCCTATACATCCGCCGTATCGTCTATCGTCGTAAACGTTGCGTAATCCGCCGGTAAACGGTAAACGTTGCGTAACCTATACATCCGCCGTCCGCCGGTAAACGTCCGCCGTATCGTCGTAAACGGTAAACGTTGCGTAATATCGTCTCCGCCGTATCGTCGTAAACGTATCGTCCCTATACAGTAAACGTTGCGTAACCTATACACCTATACATATCGTCTCCGCCGCCTATACATATCGTCTTGCGTAA", 14))
# solveFromFile()
# clumps = findClumps("CCAGCAAACACACGGCAGCAAACACAAACACACGGTTTCTTCGAAATCCCGTCTGTTAACGGGGCTGTTCGTTAACTCTGACCGTGAGGAACCCCAATCGACTGCTCTCAGTCGCAGGTACATTCTAGTCATTCGTCCGGAGGTCAAGACGGCGCCTTCGCTAATGATGATGTGGATAGCCGAGTCAAAATTGTTCGACCCCCCCCTAGCTATGTTCTTGACATTGGTTAACACAAGGTCAGTACATAAACTAGGACAGCTTAGTTTGTCTAAGCTGCAGTTACGCCAGGTTGTAAGAGTGATAGAGCTTCGAAATCGGTAGCATGCCGCCATGCGCCGGCGCATGAAGTGAGATGTTCTTTAGGATGCTCCCTGGCACACACTACTCCCACTACTAACTACTCCCCATTAACCCATTCCATTAACGAACGACCCACCCCCCAACCCCCCCGAACCCGAACCGTGTAAGAGTGTATGTCGTATTTGAGCGGCCAATGGAGCTATTTGACTATTTGACTTTATTTGACTGTAGGTGCCGTGCCGGGGCGGGGGCCGGGGTATGGGCACCGGCGGTAGCCAAAACGGTCGGTAGCCAATGTTAATCCTTCTTTGCGTGTGCCGCCCACCCCAGATAACTGGCCATAGAGGACTTATTAAGCTGGGATAAGAGAATAGCCGAGTGGTTGGAGTGGTTGGAGGTGGTTGGATCTACAGTTAGGTCCACATCCCAGCGTCCATGTTATGGATTCCGGACGTGATCTCGCCCCAGATTGTTCTCGGGTCCGACCGGCGGCGGTCGTGTCGGAGCCAGGGCAGCCTTGCTGCTCCTGCATGTAAGAAAGAAGAGCATGTAAGATCATATTTTTCTTATATTCAGGCAGCTGGGCCCTACATAGCTGTTTTATCCATAGGCAAACGCATTCGTTGCAAATCGGGTGGAGGTCGATTGCCGTCTGTTCTGCATGGTATTCTCGTTTCAAAGGAGCGGGAATAGTGCAAGCTGCTAGTGTTCATTATGCGCGCTAGGTATGTGGGCTAGTCAGGAGGCTATTTCTGGACATGCTCGTAACGCGGTGCATAATGTGTCACGCACAGGCCCCCACTAACCATGTTTCCGCCACCGTAGGTACTGCTGATGACGAGAATCAAACTGTCGGGTCCCTTTAGACGATGACAAGAGCGGGGCTGGCTCTGGCTGTGGAGCATCCTCCGAAGCGAAAAACTAACAGGATCGGCTGGGTGCTGCTGCACCGTTGCGCTTACCTGTTTTTATTACCTGTTGTTTTGCAGACCACAAAACGTTAGGCCGCGACTCGACTCGACTGGGCTAGTCGGGTCAAAGACAACTGTCTGACCATGATAAAGTATGATGTTGCGGTCGGACGTGTCGTGTGTCGGCATCACCGCGCATCTCGGCGCCTTCCGTGTTAGGACACGTGTTAATCTCACGAGGCGGATCGATTCCGGAGATACACCAAATTATTATCGCAAGTATCGCAAGTATCGCAAGCCTCGTTGAGGGGTTCTTAATATTTGCCTCTTAGGAGACTCGATGCCTGCCCTGGTAGCGTACTGCTCGTAGAATAGGATGTGCGTCCTTAGTTCCGTTACATGCGAGAGTGTACGTCAATATGTCCAAGGCGTTTCTTATTTCTTATTTTGTACTCTTATTTTGTTCGAAAGGGAGGGTAACTGCATGACCATTCCAGGCGGTTACTGAATTAACTACAGCAACTACAGCAACTACAGCAACTACAAACTACAGCGCCCAAAACAGCCAAAACAGCCAAAACAGCCAAAACAGCCAAAACAGCCAAAACAG", 9, 28, 3)
# for i in clumps:
#     print(i, end=" ")