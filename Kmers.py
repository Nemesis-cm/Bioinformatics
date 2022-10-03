TargetFile = open('ecoli.fa')
readTF = TargetFile.read()
Tseq = "".join(readTF.split())

def cMers(Tseq,k):
    kFreq = {}

    for i in range(0, len(Tseq)- k + 1):
        kmer = Tseq [i:i + k]

    if kmer in kFreq:
        kFreq[kmer] += 1
    else:   
        kFreq[kmer] = 1

    return kFreq

rf = cMers(Tseq,6)

ListRF = ([rf[sequence], sequence] for sequence in rf)

ListRF.sort(reverse = True)

print(ListRF)
