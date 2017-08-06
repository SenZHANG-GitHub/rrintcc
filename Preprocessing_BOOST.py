from itertools import islice
freq_stat = open(r"ESRD2_HKDRGWA_6445CC_Clean_Ch1-22_frq_stat.frq",mode = 'r')
ped_file = open(r"ESRD2_HKDRGWA_6445CC_Clean_Ch1-22.ped",mode = 'r')
result_file = open(r"ESRD2_BOOST.txt",mode = 'w+')
SNPs = []
line_number = -1
freq_stat_lines = freq_stat.readlines()
for line in freq_stat_lines:
    if line_number==-1:
        line_number=line_number+1
        continue
    statistics = line.split()
    A1 = statistics[2]
    A2 = statistics[3]
    SNPs.append([A1+A1,A1+A2,A2+A1,A2+A2])

index=0

while True:
    next_n_lines = list(islice(ped_file,20))
    if not next_n_lines:
        break
    res_n_lines=[]    
    for line in next_n_lines:
        if len(line)==0:
            continue
        str_list = []
        sample = line.split()

        # if not ((sample[4]=="1" or sample[4]=="2") and (sample[5]=="1" or sample[5]=="2")):
        #     continue

        if not (sample[5]=="1" or sample[5]=="2"):
            continue

        # sample[5] -> Phenotype (1: unaffected; 2: affected; -9/0: missing)
        if sample[5] == "1":
            str_list.append("0")
            
        elif sample[5] == "2":
            str_list.append("1")
        else:
            continue
            
        # # sample[4] -> Sex (1: male; 2: female; other: missing)
        # if sample[4] =="1":
        #     str_list.append("1")
        #
        # elif sample[4]=="2":
        #     str_list.append("0")
        #
        # else:
        #     continue

        # For ped file, the SNP info start  from the seven column
        j=6

        for i in range(len(SNPs)):
            if (sample[j]+sample[j+1])==SNPs[i][0]:
                str_list.append("2")
            elif (sample[j]+sample[j+1])==SNPs[i][1] or (sample[j]+sample[j+1])==SNPs[i][2]:
                str_list.append("1")
            else:
                str_list.append("0")
            
            j=j+2
        
        res_n_lines.append(' '.join(str_list))
        if index%100==0:
            print(str(index)+"\n")
        index= index+1   
    for line in res_n_lines:
        result_file.write(line+'\n')
    
freq_stat.close()
ped_file.close()
result_file.close()    
