import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
from sequenceObj import Sequence
from sklearn.metrics.cluster import adjusted_rand_score
import random
import time

#add in likelihood for the ideal situation (initialize to the true)

#cannot do ARI with different cluster lengths
#when I did ARI, got -0.11

start_time = time.time()

K = 8
# K = number of clusters
N = 91
# N = length of each sequence
S = 1415
# S = number of data points
bestOp = 10
#bestOp is the number of runs to try and find the peak
masterStoredClustering = []
#mSC stores all the clusterings
sumTestArrBest=[]
#sumTestArrBest stores all the WPSums
sortedSeqObjArr = []
finalClusters = [0]*K
# A = 0
# C = 1
# T = 2
# G = 3

# data = [[A,A],[C,C],[T,T],[G,G],[A,A],[C,C],[T,T]]
# print(data)
# data_full = [[G,A]]
# size = 5
# for i in range(size):
# data_full = np.concatenate((data_full, data))

# print(data_full)
#this code segments reads ATCG values from a flie, and converts them to arrays of 1,2,3,4
#yep = open("/Users/mallika/PycharmProjects/DirichletBio/venv/lib/shuffledExtendedseq.txt", "r")
#writtenText is unshuffled
#yep = open("/Users/mallika/PycharmProjects/DirichletBio/venv/lib/leeluMtest.txt", "r")
yep = open("/Users/mallika/PycharmProjects/DirichletBio/venv/lib/newflyData.txt", "r")
#yep = open("/Users/mallika/PycharmProjects/DirichletBio/venv/lib/mAlphashuffle.txt", "r")
finallist = []
seqlist = []
seqOBJArr = []
# seqlist is the list of all dna converted to a an array of 2 elements
for line in yep:
    values = line.split()
    seq = []
    #line = yep.readline()
    for j in range(N):
        if values[0][j] == "A":
            seq.append(0)
        elif values[0][j] == "C":
            seq.append(1)
        elif values[0][j] == "T":
            seq.append(2)
        else:
            seq.append(3)
    seqlist.append(seq)
    #s = Sequence.make_sequence(int(values[1]), seq)
    s = Sequence.make_sequence(-1, seq)
    seqOBJArr.append(s)
yep.close()
seqarray=np.asarray(seqlist)

preCluster = []
for i in range(seqOBJArr.__len__()):
    preCluster.append(seqOBJArr[i].value)



fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, figsize=(4,10))

#ax1.imshow(sortedList, extent=[0,100,0,1])
#ax1.set_title('Default')

cmap = mpl.colors.ListedColormap(['green','blue', 'red', 'orange'])



ax2.imshow(preCluster, cmap=cmap,extent=[0,100,0,1], aspect='auto')

#ax2.set_title('Fly Data Clustering')

#ax3.imshow(sortedList, extent=[0,100,0,1], aspect=100)
#ax3.set_title('Manually Set Aspect')

plt.tight_layout()
plt.show()



#looking at the trueWpSum

#first, init to the truth
'''
for i in range(seqOBJArr.__len__()):
    seqOBJArr[i].currentCluster=seqOBJArr[i].trueCluster

counter = [[]]
counter = [[1 for i in range(N)] for i in range(4*K)]

def counterCalc (seqOBJArr, counter):
    #iterating through the sequence array and entering the number of ACTG for each cluster. For example, in a scenario with four clusters, the table would have 4 times 4 = 16 rows and 2 columns (N)
    #resets all the rows and columns of the counter to one
    for row in range(4*K):
        for column in range (N):
            counter[row][column]= 1
    #recalculates the counts of each ACTG of each cluster
    for i in range (S):
        k = seqOBJArr[i].currentCluster
        for j in range (N):
            m = seqOBJArr[i].value[j]
            counter [4*k + m][j] = counter [4*k + m][j] + 1

probability = [[]]
probability =  [[0 for i in range(N)] for i in range(4*K)]
def probabilityCalc (probability, counter):
    #goes 1 column, 1 cluster at a time and uses the counter array and a sum function calculate the probability of each ACTG
        for k in range (K):
            for n in range (N):
                summ = 0
                for b in range (4):
                    summ = summ + counter [4*k + b][n]
                for b in range (4):
                    probability[4*k+b][n] = counter [4*k + b][n]/summ

counterCalc(seqOBJArr, counter)
probabilityCalc(probability,counter)
print(counter[0][0], counter[1][0], counter[2][0], counter[3][0])
print(probability[0][0], probability[1][0], probability[2][0], probability[3][0])
trialTrueFreq = [4]*K
#the following loop calculates the trial frequencies of each cluser
for a in range (seqOBJArr.__len__()):
    trialTrueFreq[seqOBJArr[a].currentCluster]+= 1
#print("trialNum", trialFreq)

trialTrueFreqsum = 0
for a in range(trialTrueFreq.__len__()):
    trialTrueFreqsum += trialTrueFreq[a]
#find sum of trial frequencies

trialTrueFreqProb=[0]*K
#make trial Freq actuall a frequency by dividing it by sum
for a in range(trialTrueFreq.__len__()):
    trialTrueFreqProb[a]= trialTrueFreq[a]/trialTrueFreqsum
print(trialTrueFreqProb)

TrueTrialProb=[[0 for i in range(K)]for j in range(S)]
for s in range (S):
    #print("s", s)
    #print("smallLen", smallSeqOBJ.__len__())
    for k in range (K):
        probcurrent=0
        for n in range (N):
            seq= seqOBJArr[s].value[n]
            p = probability [4*k+seq][n]
            probcurrent = probcurrent + math.log10(p)
        TrueTrialProb[s][k]= probcurrent

WPTruesum = 0
for l in range (TrueTrialProb.__len__()):
    for k in range (K):
        WPTruesum += (10**(TrueTrialProb [l][k])) * trialTrueFreqProb[k]
        #print("prob= ", (TrainingTrialProb [l][k]), "*", "trialFreq",trialTrainFreqProb[k], "final product", wProb, "wPsum", WPsum)
        #10 raised to to convert trial freq to probability (it is currently logs)
print("trueWPSUm", WPTruesum)
'''

for bOp in range(bestOp):
#Method B: randomly assigns each seq to a cluster
    for i in range (S):
        c = random.randint(0,(K-1))
        #assign current cluster to the objects
        seqOBJArr[i].currentCluster = c

    counter = [[]]
    counter = [[1 for i in range(N)] for i in range(4*K)]

    probability = [[]]
    probability =  [[0 for i in range(N)] for i in range(4*K)]

    def counterCalc (seqOBJArr, counter):
    #iterating through the sequence array and entering the number of ACTG for each cluster. For example, in a scenario with four clusters, the table would have 4 times 4 = 16 rows and 2 columns (N)
        #resets all the rows and columns of the counter to one
        for row in range(4*K):
            for column in range (N):
                counter[row][column]= 1
        #recalculates the counts of each ACTG of each cluster
        for i in range (S):
            k = seqOBJArr[i].currentCluster
            for j in range (N):
                m = seqOBJArr[i].value[j]
                counter [4*k + m][j] = counter [4*k + m][j] + 1

    def probabilityCalc (probability, counter):
    #goes 1 column, 1 cluster at a time and uses the counter array and a sum function calculate the probability of each ACTG
        for k in range (K):
            for n in range (N):
                summ = 0
                for b in range (4):
                    summ = summ + counter [4*k + b][n]
                for b in range (4):
                    probability[4*k+b][n] = counter [4*k + b][n]/summ


    def ClusterCalc(seqOBJArr, probability, counter):
        changeCount = 10
        while (changeCount != 0):
            changeCount = 0
            #counts the number of times a sequence's cluster is reassigned
            counterCalc(seqOBJArr, counter)
            probabilityCalc (probability, counter)
            # print("cluster: ",cluster_id)
            # print("probability: ", probability)
            for s in range (S):
                prob = -1000000000000000
                optimclust = -1
                for k in range (K):
                    probcurrent= 0
                    for n in range (N):
                        seq = seqarray[s][n]
                        p = probability [4*k+seq][n]
                        probcurrent = probcurrent + math.log10(p)
                    if (probcurrent > prob):
                        #print("probcurrent: ",probcurrent," prob: ", prob)
                        prob = probcurrent
                        optimclust = k

                old_cluster_id1 = seqOBJArr[s].currentCluster
                #print("old :", old_cluster_id1)
                seqOBJArr[s].currentCluster= optimclust
                #print("optim :", optimclust)
                #old_cluster_id = cluster_id[s]
                #print("old2 :", old_cluster_id1)
                #cluster_id[s]= optimclust
                #print("optim2 :", optimclust)
                #print(cluster_id[s])
                #print(optimclust, old_cluster_id1)
                if (optimclust !=  old_cluster_id1):
                    changeCount = changeCount + 1
                    #print("change:", changeCount)
            print("changeCount: ", changeCount)
            #checking if there were no swaps made:

    ClusterCalc(seqOBJArr, probability, counter)
    assignedClus=[]
    for i in range(seqOBJArr.__len__()):
        assignedClus.append(seqOBJArr[i].currentCluster)
    masterStoredClustering.append(assignedClus)

    #calculating the likelihood of seq existing within given starting point
    #first, calc, the trialfrequenceies of clusters

     #empty array of length K (cluster number)
    trialTrainFreq = [4]*K
    #the following loop calculates the trial frequencies of each cluser
    for a in range (seqOBJArr.__len__()):
        trialTrainFreq[seqOBJArr[a].currentCluster]+= 1
    #print("trialNum", trialFreq)

    trialTrainFreqsum = 0
    for a in range(trialTrainFreq.__len__()):
        trialTrainFreqsum += trialTrainFreq[a]
    #find sum of trial frequencies

    trialTrainFreqProb=[0]*K
    #make trial Freq actuall a frequency by dividing it by sum
    for a in range(trialTrainFreq.__len__()):
        trialTrainFreqProb[a]= trialTrainFreq[a]/trialTrainFreqsum

    TrainingTrialProb=[[0 for i in range(K)]for j in range(S)]
    for s in range (S):
        #print("s", s)
        #print("smallLen", smallSeqOBJ.__len__())
        prob = -1000000000000000
        optimclust = -1
        for k in range (K):
            probcurrent=0
            for n in range (N):
                seq= seqOBJArr[s].value[n]
                p = probability [4*k+seq][n]
                probcurrent+= math.log10(p)
            TrainingTrialProb[s][k]= probcurrent

    trainingWP = [0]*S
    for l in range (TrainingTrialProb.__len__()):
        wProb = 0
        for k in range (K):
            #print(10**(TrainingTrialProb [l][k]))
            wProb += (10**(TrainingTrialProb [l][k])) * trialTrainFreqProb[k]
            #print("prob= ", (TrainingTrialProb [l][k]), "*", "trialFreq",trialTrainFreqProb[k], "final product", wProb, "wPsum", WPsum)
             #10 raised to to convert trial freq to probability (it is currently logs)
        trainingWP[l]= math.log10(wProb)
    logProduct = 0
    for i in range(trainingWP.__len__()):
        logProduct+= trainingWP[i]
    sumTestArrBest.append(logProduct)



#finds the greatest of the WP Sum, to find the best clustering
greatest = -100000000000000000000000000000000000
greatestPos = -1
for i in range(sumTestArrBest.__len__()):
    if (sumTestArrBest[i]> greatest):
        greatest= sumTestArrBest[i]
        greatestPos=i

bestAssignedClus= masterStoredClustering[greatestPos]


#ularValue = [0,3,5,6,1,2,4,7]

def orderByCluster(bAC, seqobjarr):
    sortedValues=[]
    for k in range (K):
    #for i in range (particularValue.__len__()):
        #k= particularValue[i]
        bool = True
        for i in range (bAC.__len__()):
            if (bool == True):
                sortedValues.append([5]*N)
                sortedValues.append([5]*N)
                sortedValues.append([5]*N)
                sortedValues.append([5]*N)
                sortedValues.append([5]*N)
                sortedValues.append([5]*N)
                sortedValues.append([5]*N)
                sortedValues.append([5]*N)
                sortedValues.append([5]*N)
                sortedValues.append([5]*N)
                sortedValues.append([5]*N)
                bool = False
            if (bAC[i] == k):
                sortedSeqObjArr.append(seqobjarr[i])
                sortedValues.append(seqobjarr[i].value)
                finalClusters[k]+=1
    return(sortedValues)

sortedVal= orderByCluster(bestAssignedClus, seqOBJArr)

#print(sortedValues)




fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, figsize=(4,10))

#ax1.imshow(sortedList, extent=[0,100,0,1])
#ax1.set_title('Default')

cmap = mpl.colors.ListedColormap(['green','blue', 'red', 'orange', 'black'])



ax2.imshow(sortedVal, cmap=cmap,extent=[0,100,0,1], aspect='auto')

#ax2.set_title('Fly Data Clustering')

#ax3.imshow(sortedList, extent=[0,100,0,1], aspect=100)
#ax3.set_title('Manually Set Aspect')

plt.tight_layout()
plt.show()
print("finalCluster", finalClusters)


trueClustersSize = [0]*10
for i in range(seqOBJArr.__len__()):
    trueClustersSize[seqOBJArr[i].trueCluster]+=1
print("trueClusters", trueClustersSize)

print("greatestPos", greatestPos)

trueClusters = []
for i in range(seqOBJArr.__len__()):
    trueClusters.append(seqOBJArr[i].trueCluster)


print("ARI", adjusted_rand_score(trueClusters, bestAssignedClus))


print("--- %s seconds ---" % (time.time() - start_time))
