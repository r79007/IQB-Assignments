alpha_dict = {'E': 1.53, 'A': 1.45, 'L': 1.34, 'H': 1.24, 'M': 1.20, 'Q': 1.17, 'W': 1.14, 'V': 1.14, 'F': 1.12, 'K': 1.07,
         'I': 1.00, 'D': 0.98, 'T': 0.82, 'S': 0.79, 'R': 0.79, 'C': 0.77, 'N': 0.73, 'Y': 0.61, 'P': 0.59, 'G': 0.53}
beta_dict = {'M': 1.67, 'V': 1.65, 'I': 1.60, 'C': 1.30, 'Y': 1.29, 'F': 1.28, 'Q': 1.23, 'L': 1.22, 'T': 1.20, 'W': 1.19,
         'A': 0.97, 'R': 0.90, 'G': 0.81, 'D': 0.80, 'K': 0.74, 'S': 0.72, 'H': 0.71, 'N': 0.65, 'P': 0.62, 'E': 0.26}

protein_seq = "SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDTVYCPRHVICTAEDMLNPNYEDLLIRKSNHSFLVQAGNVQLRVIGHSMQNCLLRLKVDTSNPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNHTIKGSFLNGSCGSVGF"

helix=[]
strand=[]

lst=[]

lst.append(0)
lst.append(1)
lst.append(6)

for i in range(len(protein_seq)):
    helix.append("#")
    strand.append("#")


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

        k=i
        while(k>=0):
            score=0
            for j in range(0,4):
                residue=alpha_dict[k-j]
                score+=alpha_dict[residue]

            if(score>=4):
                for j in range(0,4):
                    helix[k-j]="H"
                k-=1
            else:
                break











