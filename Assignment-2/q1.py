alpha_dict = {'E': 1.53, 'A': 1.45, 'L': 1.34, 'H': 1.24, 'M': 1.20, 'Q': 1.17, 'W': 1.14, 'V': 1.14, 'F': 1.12, 'K': 1.07,
         'I': 1.00, 'D': 0.98, 'T': 0.82, 'S': 0.79, 'R': 0.79, 'C': 0.77, 'N': 0.73, 'Y': 0.61, 'P': 0.59, 'G': 0.53}
beta_dict = {'M': 1.67, 'V': 1.65, 'I': 1.60, 'C': 1.30, 'Y': 1.29, 'F': 1.28, 'Q': 1.23, 'L': 1.22, 'T': 1.20, 'W': 1.19,
         'A': 0.97, 'R': 0.90, 'G': 0.81, 'D': 0.80, 'K': 0.74, 'S': 0.72, 'H': 0.71, 'N': 0.65, 'P': 0.62, 'E': 0.26}

protein_seq = "SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDTVYCPRHVICTAEDMLNPNYEDLLIRKSNHSFLVQAGNVQLRVIGHSMQNCLLRLKVDTSNPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNHTIKGSFLNGSCGSVGF"

helix=[]
strand=[]

result=[]

lst=[]

lst.append(0)
lst.append(1)
lst.append(6)

for i in range(len(protein_seq)):
    helix.append("#")
    strand.append("#")
    result.append("#")

#For helix

i=0
while(i<len(protein_seq)-5):
    #Counting Favourable residues
    num_fav=lst[0]
    for j in range(lst[2]):
        residue=protein_seq[i+j]
        if(alpha_dict[residue]>=1):
            num_fav+=lst[1]

    #

    if(num_fav>=4):
        for j in range(lst[2]):
            helix[i+j]="H"

        #Now we do extensions

        k=i+3
        while(k-3>=0):
            score=0
            for j in range(0,4):
                residue=protein_seq[k-j]
                score+=alpha_dict[residue]

            if(score>=4):
                for j in range(0,4):
                    helix[k-j]="H"
                k-=1
            else:
                break

        k = i+2
        while (k+4 <= len(protein_seq)):
            score = 0
            for j in range(0, 4):
                residue = protein_seq[k + j]
                score += alpha_dict[residue]

            if (score >= 4):
                for j in range(0, 4):
                    helix[k + j] = "H"
                k += 1
            else:
                break

    i+=1




#For strand

i=0

while(i<len(protein_seq)-4):

    # Counting Favourable residues
    num_fav = lst[0]
    for j in range(lst[2]-1):
        residue = protein_seq[i + j]
        if (beta_dict[residue] >= 1):
            num_fav += lst[1]

    #

    if (num_fav >= 3):
        for j in range(lst[2]-1):
            strand[i + j] = "S"

        # Now we do extensions

        k = i + 2
        while (k - 3 >= 0):
            score = 0
            for j in range(0, 4):
                residue = protein_seq[k - j]
                score += beta_dict[residue]

            if (score > 4):
                for j in range(0, 4):
                    strand[k - j] = "S"
                k -= 1
            else:
                break

        k = i + 2
        while (k + 4 <= len(protein_seq)):
            score = 0
            for j in range(0, 4):
                residue = protein_seq[k + j]
                score += beta_dict[residue]

            if (score > 4):
                for j in range(0, 4):
                    strand[k + j] = "S"
                k += 1
            else:
                break

    i+=1

i=0



while(i<len(protein_seq)):

    idx=0

    if(helix[i]=='#' and strand[i]=='#'):
        result[i]='T'
        i+=1

    elif(helix[i]=='H' and strand[i]=='#'):
        result[i]='H'
        i+=1

    elif(helix[i]=='#' and strand[i]=='S'):
        result[i]='S'
        i+=1

    else:

        while(idx+i<len(protein_seq)):
            if(helix[i+idx]=='H' and strand[i+idx]=='S'):
                idx+=1
            else:
                break

        hel_score = 0
        strand_score = 0

        for j in range(0, idx):
            residue_alpha=protein_seq[i+j]
            residue_beta = protein_seq[i + j]
            hel_score+=alpha_dict[residue_alpha]
            strand_score+=beta_dict[residue_beta]

        if(strand_score>hel_score):
            for j in range(0, idx):
                result[i+j]='S'
        else:
            for j in range(0, idx):
                result[i + j] = 'H'

        i+=idx

server_result="TTTT     HHHHHH EEEEEETTEEEEEEEETTEEEEEGGGG  HHHHH   HHHHHHH  GGG EEEETTEEE EEEEEEETTEEEEEE   TTTT        TTTEEEEEEEEETTEEEEEEEEEETTTT B    TTTTTTTEE "
result=''.join(result)
print(protein_seq)
print(result)
print(server_result)

print()

print("Differing Regions(The first string represents the server result, and the second string represents my result): ")
print("Here the gaps in the server result have been substituted by C, because the gaps actually mean coils")

idx=0

change = False

diff_server_result=""
diff_result=""

while(idx<len(protein_seq)):
    i = idx
    change = False
    while(i<len(protein_seq) and (server_result[i]!=result[i])):
        if(not(server_result[i]=="E" and result[i]=="S")):

            if(result[i]==" "):
                diff_result+="C"
            else:
                diff_result += result[i]
            if(server_result[i]==" "):
                diff_server_result+="C"
            else:
                diff_server_result += server_result[i]




            # print(server_result[i], end="")
            # print(result[i], end="")
            change=True
            i+=1


        else:
            break




    if(change==True):

        print()

        print("Indices: ", idx, " to ", i - 1)

        print()

        print(diff_server_result[idx-1:i])
        print(diff_result[idx-1:i])



        idx=i
    else:
        idx += 1
        diff_result+=" "
        diff_server_result+=" "

# print(diff_server_result)
# print(diff_result)





























