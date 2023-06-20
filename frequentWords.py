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

print(frequentWords("ACGTTGCATGTCGCATGATGCATGAGAGCT", 3))
print(frequencyTable("ACGTTGCATGTCGCATGATGCATGAGAGCT", 3))
print(betterFrequentWords("CCTATACACCTATACATTGCGTAACCTATACACCTATACATCCGCCGTATCGTCTATCGTCTTGCGTAACCTATACATCCGCCGTCCGCCGGTAAACGTCCGCCGTCCGCCGTCCGCCGTATCGTCCCTATACATATCGTCTTGCGTAATCCGCCGCCTATACATATCGTCGTAAACGTCCGCCGGTAAACGTTGCGTAACCTATACAGTAAACGTATCGTCTCCGCCGTATCGTCTTGCGTAATATCGTCCCTATACAGTAAACGTATCGTCTCCGCCGCCTATACATATCGTCCCTATACAGTAAACGGTAAACGTATCGTCTCCGCCGTTGCGTAACCTATACATATCGTCTTGCGTAATATCGTCCCTATACATTGCGTAATCCGCCGTATCGTCTTGCGTAACCTATACATATCGTCCCTATACAGTAAACGGTAAACGTTGCGTAACCTATACATATCGTCGTAAACGTCCGCCGTCCGCCGCCTATACACCTATACATCCGCCGCCTATACAGTAAACGGTAAACGCCTATACATATCGTCTCCGCCGCCTATACATATCGTCTCCGCCGCCTATACATATCGTCTTGCGTAATATCGTCCCTATACATCCGCCGTATCGTCTATCGTCGTAAACGTTGCGTAATCCGCCGGTAAACGGTAAACGTTGCGTAACCTATACATCCGCCGTCCGCCGGTAAACGTCCGCCGTATCGTCGTAAACGGTAAACGTTGCGTAATATCGTCTCCGCCGTATCGTCGTAAACGTATCGTCCCTATACAGTAAACGTTGCGTAACCTATACACCTATACATATCGTCTCCGCCGCCTATACATATCGTCTTGCGTAA", 14))