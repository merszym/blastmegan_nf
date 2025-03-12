from pysam import AlignmentFile
from Bio.Seq import Seq

with open("pseudouniq/pseudouniq_stats.txt",'a') as stats:
    for readfile,writefile in zip (input, output):
        samfile = AlignmentFile(readfile, "rb")
        
        #using template=samfile to copy the header
        bamout = AlignmentFile(writefile, "wb", template=samfile)
        #open in appending mode
        #stats = open("pseudouniq/pseudouniq_stats.txt",'a')
        duplicate = {}
        length_passed = 0
        passed = 0
        allreads = 0
        
        for aln in samfile.fetch(until_eof=True):
            allreads += 1
            # flag 5 is bit 1 and 4 (-F up)
            if((not (aln.flag & 5)) and (aln.query_length >= rm_dup_query_len)): #35
                length_passed += 1
                if aln.is_reverse: # if the read is reverse
                    # the sequence has to be revcomp'ed'
                    seq = Seq(aln.query_sequence).reverse_complement()#create a sequence object
                else:
                    seq = aln.query_sequence
                try: #if key doesn't exists, create in except below
                    #duplicate[aln.query_alignment_sequence] += 1
                    duplicate[seq] += 1
                    #if(aln.query_sequence == "CGTGCGTACACGTGCGTACACGTGCGTACACGTGCG"):
                    #    print(duplicate[aln.query_alignment_sequence],aln)
                except KeyError:
                    duplicate[seq] = 1
                    #duplicate[aln.query_alignment_sequence] = 1
                    #if(aln.query_sequence == "CGTGCGTACACGTGCGTACACGTGCGTACACGTGCG"):
                    #    print(aln)
                #keep only reads apearing at least twice
                #twice by defaul, now can be assigned through command line
                # we want the sequence to be printed only once it reaches the min
                #number of duplicates
                if(duplicate[seq] == num_min_dup):
                    passed += 1
                    bamout.write(aln)

        samfile.close()
        bamout.close()        
        #write stats
        dups_average = [v for v in duplicate.values() if (v>=num_min_dup)]
        all_average = [v for v in duplicate.values()]
        #avoid division by 0
        if (len(all_average)):
            rounded_avrg_all = "{:.1f}".format(sum(all_average)/len(all_average))
        else:
            rounded_avrg_all = "0"
        if(len(dups_average)):
            rounded_avrg_final =" {:.1f}".format(sum(dups_average)/len(dups_average))
        else:
            rounded_avrg_final = "0"
        print ("#file","all_seqs","seqs>="+str(rm_dup_query_len), \
        "avrg_times_seen_L>="+str(rm_dup_query_len), \
            "final_noPCRdups_seqs","avrg_times_seen_final",sep="\t",file=stats) #header
        print (wildcards["sample"], allreads, length_passed, rounded_avrg_all,\
            passed, rounded_avrg_final,sep="\t",file=stats) #print out the stats