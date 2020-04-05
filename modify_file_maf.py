import sys

file_2_modify=open(sys.argv[1],'r')
file_2_write=open(sys.argv[2],'w')
species1=sys.argv[3]
species2=sys.argv[4]
for line in file_2_modify:
    elem=line.split('\t')
    file_2_write.write(elem[0]+'\t'+str(species1+'.'+elem[1])+'\t'+elem[2]+'\t'+str(species2+'.'+elem[3])+'\t'+elem[4]+'\t'+elem[5]+'\t'+elem[6])

file_2_modify.close()
file_2_write.close()