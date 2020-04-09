import sys

file_2_modify=open(sys.argv[1],'r')
file_2_write=open(sys.argv[2],'w')
count=0
for line in file_2_modify:
    count+=1
    file_2_write.write(str(count)+'\t'+line)

file_2_modify.close()
file_2_write.close()