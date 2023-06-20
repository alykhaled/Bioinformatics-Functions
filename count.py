def patternCount(text, pattern):
    """
    Counts the number of times a pattern appears in a text
    :param text: the text to search
    :param pattern: the pattern to search for
    :return: the number of times the pattern appears in the text

    For example to count the number of times "ATAT" appears in "GATATATGCATATACTT"
    """
    count = 0
    for i in range(len(text) - len(pattern) + 1):
        if text[i:i+len(pattern)] == pattern:
            count += 1
    return count

print(patternCount("CTGCGCCCGTCCCCCGTCCTTCCCGTCCCCCGTCCATTCGCCCGTCCGGCCCGTCCCCCCCGTCCCCCGTCCATCCCGTCCTAGGGCCCGTCCTCCGACC", "CCCGTCCCC"))
    