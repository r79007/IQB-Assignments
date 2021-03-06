alpha_dict = {'E': 1.53, 'A': 1.45, 'L': 1.34, 'H': 1.24, 'M': 1.20, 'Q': 1.17, 'W': 1.14, 'V': 1.14, 'F': 1.12, 'K': 1.07,'I': 1.00, 'D': 0.98, 'T': 0.82, 'S': 0.79, 'R': 0.79, 'C': 0.77, 'N': 0.73, 'Y': 0.61, 'P': 0.59, 'G': 0.53}
beta_dict = {'M': 1.67, 'V': 1.65, 'I': 1.60, 'C': 1.30, 'Y': 1.29, 'F': 1.28, 'Q': 1.23, 'L': 1.22, 'T': 1.20, 'W': 1.19,'A': 0.97, 'R': 0.90, 'G': 0.81, 'D': 0.80, 'K': 0.74, 'S': 0.72, 'H': 0.71, 'N': 0.65, 'P': 0.62, 'E': 0.26}

protein_seq = "SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDTVYCPRHVICTAEDMLNPNYEDLLIRKSNHSFLVQAGNVQLRVIGHSMQNCLLRLKVDTSNPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNHTIKGSFLNGSCGSVGF"

#Helix list represents the sites predicted as part of Helix
#Strand list represents the sites predicted as part of Strand
#Result list will contain our final answer

helix=[]
strand=[]

result=[]

lst=[]

#lst and common_arr store the common numbers that will be used(but could not implement it in a useful way, so there is no real use of it)
common_arr=[]

common_arr.append(2)
common_arr.append(3)
common_arr.append(4)
common_arr.append(5)

lst.append(0)
lst.append(1)
lst.append(6)

for i in range(len(protein_seq)):
    helix.append("#")
    strand.append("#")
    result.append("#")
#Wherever I have used the word 'score', it means average propensity
#For helix

i=0
#"We check all possible stretches of length 6
end_idx=len(protein_seq)-5
while(i<end_idx):
    #Counting Favourable residues in these stretches
    num_fav=lst[0]
    #Checking in stretches of size 6
    for j in range(lst[2]):
        residue=protein_seq[i+j]
        #If P(H)>=1 then we increment num_fav(favourable residues)
        if(alpha_dict[residue]>=1):
            num_fav+=lst[1]

    #If then number of favourable residues is >=4 then this site can have a helix

    if(num_fav>=common_arr[2]):
        for j in range(lst[2]):
            helix[i+j]="H"

        #Now we do extensions

        #We first try to extend in stretches of 4 to the right
        k=i+common_arr[1]
        while(k>=3):
            score=0
            for j in range(0,4):
                site_idx = k - j
                residue=protein_seq[site_idx]
                score+=alpha_dict[residue]
            #If the current right stretch has sum of the P(H) of the residues >=4 then this can be a helix
            if(score>=common_arr[2]):
                for j in range(0,4):
                    site_idx = k - j
                    helix[site_idx]="H"
                k-=1

            #If not then we break
            else:
                break
        #Now we try to extend in stretches of 4 to the left
        k = i+common_arr[0]
        end_i=len(protein_seq) - 4
        while (k <= end_i):
            score = 0
            for j in range(0, 4):
                site_idx=k+j
                residue = protein_seq[site_idx]
                score += alpha_dict[residue]
            # If the current right stretch has sum of the P(H) of the residues >=4 then this can be a helix
            if (score >= common_arr[2]):
                for j in range(0, 4):
                    site_idx = k + j
                    helix[site_idx] = "H"
                k += 1
            #Else we break
            else:
                break

    i+=1




#For strand

i=0

#We check all possible stretches of length 5
end_idx=len(protein_seq)-4
while(i<end_idx):

    # Counting Favourable residues in these stretches
    num_fav = lst[0]
    # Checking in stretches of size 5
    for j in range(lst[2]-1):
        residue = protein_seq[i + j]
        # If P(H)>=1 then we increment num_fav(favourable residues)
        if (beta_dict[residue] >= 1):
            num_fav += lst[1]


    #If number of favourable residues is >= 3 then this site can be a strand
    if (num_fav >= 3):
        for j in range(lst[2]-1):
            strand[i + j] = "S"

        # Now we do extensions

        # We first try to extend in stretches of 3 to the right
        st_idx = i + common_arr[0]
        while (st_idx >= 3):
            score = 0
            for j in range(0, 4):
                site_idx=st_idx - j
                residue = protein_seq[site_idx]
                score += beta_dict[residue]
            # If sum of P(S) of the residues in this stretch is >= 4 then this site can be a strand
            if (score >= common_arr[2]):
                for j in range(0, 4):
                    site_idx=st_idx - j
                    strand[site_idx] = "S"
                st_idx -= 1
            #Else we break from working on this stretch
            else:
                break

        # Now we try to extend in stretches of 3 to the left
        st_idx = i + common_arr[0]
        end_i=len(protein_seq) - common_arr[2]
        while (st_idx  <= end_i):
            score = 0
            for j in range(0, 4):
                site_idx=st_idx + j
                residue = protein_seq[site_idx]
                score += beta_dict[residue]
            # If sum of P(S) of the residues in this stretch is >= 4 then this site can be a strand
            if (score >= common_arr[2]):
                for j in range(0, 4):
                    site_idx=st_idx + j
                    strand[site_idx] = "S"
                st_idx += 1
            # Else we break from working on this stretch
            else:
                break

    i+=1

i=0

# Now we fill in the predictions for generating the final result

while(i<len(protein_seq)):

    idx=0
    #If our method could not predict whether it could be a strand or a helix, we take that to be a blank(i.e #)
    if(helix[i]=='#' and strand[i]=='#'):
        result[i]='#'
        i+=1
    #If it is predicted as Helix but not as Strand then it could be a helix
    elif(helix[i]=='H' and strand[i]=='#'):
        result[i]='H'
        i+=1
    # If it is predicted as Strand but not as Helix then it could be a Strand
    elif(helix[i]=='#' and strand[i]=='S'):
        result[i]='S'
        i+=1
    #If it is predicted both as a strand and a helix then we check whose score is higher, the one with the higher score(sum of probabilities) has more probability to be a correct prediction
    else:

        while(idx+i<len(protein_seq)):
            #We find the limit(last index in the stretch) to which the predictions are both H and S
            if(helix[i+idx]=='H' and strand[i+idx]=='S'):
                idx+=1
            else:
                break

        hel_score = 0
        strand_score = 0

        for j in range(0, idx):
            # Now we calculate the score of helix and strand till this index
            residue_alpha=protein_seq[i+j]
            residue_beta = protein_seq[i + j]
            hel_score+=alpha_dict[residue_alpha]
            strand_score+=beta_dict[residue_beta]
        #If strand has higher score then it has more probability to be correct
        if(strand_score>hel_score):
            for j in range(0, idx):
                result[i+j]='S'
        #Else helix has more probabilty to be correct
        else:
            for j in range(0, idx):
                result[i + j] = 'H'

        i+=idx

server_result="TTTT     HHHHHH EEEEEETTEEEEEEEETTEEEEEGGGG  HHHHH   HHHHHHH  GGG EEEETTEEE EEEEEEETTEEEEEE   TTTT        TTTEEEEEEEEETTEEEEEEEEEETTTT B    TTTTTTTEE "
result=''.join(result)

print("The following sequences are the original protein sequence, the result that we got from the Chau-Fasman Method and the result from the server respectively\n")

print(protein_seq)
print(result)
print(server_result)

print()
print("Differing Regions(The first string represents the server result, and the second string represents my result): ")
print("Here the gaps in the server result have been substituted by #, they can be coils")

idx=0

change = False

diff_server_result=""
diff_result=""
#To print the differing regions
while(idx<len(protein_seq)):
    i = idx
    change = False
    while(i<len(protein_seq) and (server_result[i]!=result[i])):
        if(not(server_result[i]=="E" and result[i]=="S")):

            if(result[i]==" "):
                diff_result+="#"
            else:
                diff_result += result[i]
            if(server_result[i]==" "):
                diff_server_result+="#"
            else:
                diff_server_result += server_result[i]


            change=True
            i+=1


        else:
            break




    if(change==True):

        print()

        print("Indices: ", idx, " to ", i - 1)

        print()

        print(diff_server_result[idx:i])
        print(diff_result[idx:i])



        idx=i
    else:
        idx += 1
        diff_result+=" "
        diff_server_result+=" "

# print(diff_server_result)
# print(diff_result)





























