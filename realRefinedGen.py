import numpy as np
import random
import math

#this file generates clusters randomly, ie, creates a dirichlet distribution and creates clusters in that proportion

a = np.array([1, 2, 3])
a = np.array([1, 2, 3])
import matplotlib

matplotlib.use('TkAgg')

import matplotlib.pyplot as plt



# alpha distribution
alphaGen = 1
totalSeqNum = 1415
numClusters = 10
x = 0.2
length = 90
refArray = ["A", "C", "T", "G"]
masterCluster = []


# s = np.array([[0.250226366, 0.0002297601832, 0.165059640, 0.584416391],
# [0.0000194563325, 0.000000000245599067, 0.337103929, 0.662876614]])
# s = np.array([[2.50226366e-01, 2.97601832e-04, 1.65059640e-01, 5.84416391e-01],
# [1.94563325e-05, 2.45599067e-10, 3.37103929e-01, 6.62876614e-01]])


# printing a and b prob. of first cluster as one array

# setting a and b prob. of first cluster

#generate the distribution of clusters- basically what proportion of the total that each should be
g = np.random.dirichlet([alphaGen]*numClusters, 1)
arr = g[0]
actualNums = [0]*numClusters
total = 0
for i in range(actualNums.__len__()):
    if (i < (actualNums.__len__()-1)):
        actualNums[i]= math.floor(arr[i]*totalSeqNum)
        total += actualNums[i]

    else:
        actualNums[i]= totalSeqNum-total


print("dirichlet", arr)
print("actualNums", actualNums)

# printing aI prob
#clearing the written file and the shuffle file
f = open('/Users/mallika/PycharmProjects/DirichletBio/venv/lib/writtenFile.txt', 'r+')
f.truncate(0)
f.close()

f = open('/Users/mallika/PycharmProjects/DirichletBio/venv/lib/shuffledwrittenFile.txt', 'r+')
f.truncate(0)
f.close()


f = open("/Users/mallika/PycharmProjects/DirichletBio/venv/lib/writtenFile.txt", "w+")


for m in range(numClusters):
    s = np.random.dirichlet((x, x, x, x), length)
    allSequences = []
    for i in range(actualNums[m]):
        currentSeq = []
        for j in range(length):
            # generating the random data
            # print("")
            randA = random.random()
            # print("apoint rand is: ", randA)
            # should it be less than equal to?
            # sometimes data is printed in decimal sometimes with e-why?
            pointProb = s[j]
            if randA < pointProb[0]:
                a = 0
            elif randA < (pointProb[0] + pointProb[1]):
                a = 1
            elif randA < (pointProb[0] + pointProb[1] + pointProb[2]):
                a = 2
            elif randA <= (pointProb[0] + pointProb[1] + pointProb[2] + pointProb[3]):
                a = 3
            currentSeq.append(refArray[a])
            print(refArray[a], end="")
            #writes each letter of the sequence
            f.write(refArray[a])

        print(" ",  m, "")
        #puts space cluster and goes to anew line in the file f
        print(" ",  m, "", file = f)
        allSequences.append(currentSeq)
    masterCluster.append(allSequences)
    print("al;lseq :", allSequences.__len__())


'''
#steps to refresh simulation
- run realRefined gen
- run shuffler
- run plotting 2.0
'''

