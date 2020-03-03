#!/usr/bin/python2
#usage: python2 /filter_list_from_list.py /Candidate_snps.txt /SNP_names_from_vcf.txt /snp_codes_list.txt /batch_0.catalog.tags.tsv > /tags_candidate_snps.txt

#Get file with list of candidate SNPs (and respective POS and ID) and gets sequence tags (catalog tags from Stacks) of candidate SNPs

from sys import argv

file_1=open(argv[1],"r") # open file Candidate_snps.txt containing list of candidate snps (ID of SNPs is the order of the snp in the file .geste)

candidate_snp_list=[]
for line in file_1:
    if line.startswith("#"):
        pass
    else:
        candidate_snp=line.split("\t")[0].strip("\n") # ID of the candidate SNPs (order of the snp in the file .geste)
        candidate_snp_list.append(candidate_snp)
file_1.close()

candidate_snp_set=set(candidate_snp_list)

file_2=open(argv[2],"r") # open file SNP_names_from_vcf.txt with all SNPs and respective POS and ID from Stacks

write_file=open(argv[3], "a") # create file snp_codes_list.txt

tag_id_list=[]

for line in file_2:
    if line.startswith("#"):
        pass
    else:
        snp_order=line.split("\t")[0] # ID of the candidate SNPs (order of the snp in the file .geste)
        snp_pos=line.split("\t")[1] # snp POS (from Stacks)
        tag_id=line.split("\t")[2] # tag ID (from Stacks)
        if snp_order in candidate_snp_set:
            tag_id_list.append(tag_id)
            write_file.write(snp_order + "\t" + snp_pos +"\t" +tag_id+"\n")
file_2.close()
tag_id_set=set(tag_id_list)

catalog_file=open(argv[4],"r") # open catalog tag file from Stacks

for line in catalog_file:
    if line.startswith("#"):
        pass
    else:
        snp_name=line.split("\t")[2] # snp name in .geste
        sequence_tag=line.split("\t")[9] # sequence tag
        if snp_name in tag_id_set:
            print ">"+str(snp_name)
            print str(sequence_tag[11:]) # remove 11 first bases (6 for SbfI recognition site and 5 for barcode)

catalog_file.close() 
