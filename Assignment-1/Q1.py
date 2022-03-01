match=2
mismatch=-1
gap=-1

def initialize_matrix(seq1, seq2):
    score=0
    global maxScore
    maxScore=0
    
    mat= [[0 for i in range(len(seq2)+1)] for j in range(len(seq1) + 1)]

    for i in range (0,len(seq1)+1):
        mat[i][0]=score
        score+=gap
    
    score=0
    for i in range (0,len(seq2)+1):
        mat[0][i]=score
        score+=gap



    for i in range(1, len(seq1)+1):
        for j in range(1,len(seq2)+1):
            if(seq1[i-1]==seq2[j-1]):
                mat[i][j]=match+mat[i-1][j-1]
                if(mat[i][j]>maxScore):
                    maxScore=mat[i][j]
            else:
                mat[i][j]=max(mat[i-1][j]+gap, mat[i][j-1]+gap, mat[i-1][j-1]+mismatch)
                if(mat[i][j]>maxScore):
                    maxScore=mat[i][j]
            

    return mat


def print_matrix(A,seq1,seq2):

    for i in range (0,len(seq1)+1):
        for j in range(0, len(seq2)+1):
            print(A[i][j], end=" ")

        print()

opt1=[]
opt2=[]

def findAllOptAlign(seq1, seq2, i, j, newSeq1, newSeq2, currScore, arr):

    if(currScore==maxScore and (i==0 and j==0) or (i<0 or j<0)):
        opt1.append(newSeq1)
        opt2.append(newSeq2)
        return

    if((i==0 and j==0) or (i<0 or j<0)):
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



def print_opt_alignments(seq1, seq2, arr):


    findAllOptAlign(seq1,seq2, len(seq1), len(seq2), "", "", 0,arr)

    idx=-1
    if(len(opt1)>1):
            
            print("There are more than 1 optimal alignments possible\n")
    else:

            print("There is only 1 optimal alignment possible\n")

    print()
    for j in range(len(opt1)):
            print(opt1[j][::-1])
            print(opt2[j][::-1])
            print('Score: ',maxScore)
            print('-------------------------------')

    

seq2 = "ATCAGAGTA"
seq1= "TTCAGTA"

arr=initialize_matrix(seq1, seq2)
print("The bi-dimensional array is \n")

print_matrix(arr,seq1,seq2)
print('\n')

print_opt_alignments(seq1, seq2, arr)







