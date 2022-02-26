
match=2
mismatch=-1
gap=-2

seq2 = "ATCAGAGTA"
seq1= "TTCAGTA"


def initialize_matrix():
    score=0
    
    mat= [[0 for i in range(len(seq2)+1)] for j in range(len(seq1) + 1)]

    for i in range (0,len(seq1)+1):
        mat[i][0]=0
        
    
    for i in range (0,len(seq2)+1):
        mat[0][i]=0

    for i in range(1, len(seq1)+1):
        for j in range(1,len(seq2)+1):
            if(seq1[i-1]==seq2[j-1]):
                mat[i][j]=match+mat[i-1][j-1]
            else:
                mat[i][j]=max(mat[i-1][j]+gap, mat[i][j-1]+gap, mat[i-1][j-1]+mismatch,0)
            

    return mat


def print_matrix(A):

    for i in range (0,len(seq1)+1):
        for j in range(0, len(seq2)+1):
            print(A[i][j], end=" ")

        print()

opt1=[]
opt2=[]
scores=[]

def findAllOptAlign(seq1, seq2, i, j, newSeq1, newSeq2, currScore, arr):

    if(arr[i][j]==0):
        opt1.append(newSeq1)
        opt2.append(newSeq2)
        scores.append(currScore)
        return
        

    if((i==0 and j==0) or (i<0 or j<0)):
        opt1.append(newSeq1)
        opt2.append(newSeq2)
        scores.append(currScore)
        return
    
    if(i==0):
        findAllOptAlign(seq1, seq2, i, j-1, newSeq1+"-", newSeq2+seq2[j-1], currScore+gap, arr)
        return
    if(j==0):

        findAllOptAlign(seq1, seq2, i-1, j, newSeq1+seq1[i-1], newSeq2+"-", currScore+gap, arr)
        return
    
    #now we try to find the optimal alignments from all three possiblities -> match, mismatch, or gap

    if(seq1[i-1]==seq2[j-1]):
        findAllOptAlign(seq1, seq2,i-1,j-1,newSeq1+seq1[i-1], newSeq2+seq2[j-1], currScore+match, arr)
        if(arr[i][j]-gap==arr[i-1][j]):
            findAllOptAlign(seq1, seq2, i-1, j, newSeq1+seq1[i-1], newSeq2+"-", currScore+gap, arr)
        if(arr[i][j]-gap==arr[i][j-1]):
            findAllOptAlign(seq1, seq2, i, j-1, newSeq1+"-", newSeq2+seq2[j-1], currScore+gap, arr)

    else:
        findAllOptAlign(seq1, seq2,i-1,j-1,newSeq1+seq1[i-1],newSeq2+seq2[j-1],currScore+mismatch, arr)
        if(arr[i][j]-gap==arr[i-1][j]):
            findAllOptAlign(seq1, seq2, i-1, j, newSeq1+seq1[i-1], newSeq2+"-", currScore+gap, arr)
        if(arr[i][j]-gap==arr[i][j-1]):
            findAllOptAlign(seq1, seq2, i, j-1, newSeq1+"-", newSeq2+seq2[j-1], currScore+gap, arr)


def find_max_indices(arr):

    max_score=float('-inf')
    max_i=-1
    max_j=-1

    for i in range(len(seq1)+1):
        for j in range(len(seq2)+1):
            if(max_score<arr[i][j]):
                max_score=arr[i][j]
                max_i=i
                max_j=j

    return (max_i,max_j)

def print_opt_alignments(seq1, seq2, arr):

    (idx_i, idx_j)=find_max_indices(arr)

    findAllOptAlign(seq1,seq2, idx_i, idx_j, "", "", 0,arr)

    idx=-1
    maxScore=float('-inf')


    for i in scores:
        if(i>maxScore):
            maxScore=i

    for j in range(len(opt1)):
        if(scores[j]==maxScore):
            print(opt1[j][::-1])
            print(opt2[j][::-1])
            print(scores[j])
            print('-------------------------------')



arr=initialize_matrix()
print("The bi-dimensional array is\n")
print_matrix(arr)
print('\n')
print("There is only 1 optimal alignment possible\n")
print_opt_alignments(seq1, seq2, arr)