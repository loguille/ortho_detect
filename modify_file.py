import sys

file_2_modify=open(sys.argv[1],'r')
file_2_write=open(sys.argv[2],'w')
species=sys.argv[3]
for line in file_2_modify:
    elem=line.split('\t')
    file_2_write.write(elem[0]+'\t'+str(species+'.'+elem[1])+'\t'+elem[2]+'\t'+elem[3]+'\t'+elem[4])

file_2_modify.close()
file_2_write.close()