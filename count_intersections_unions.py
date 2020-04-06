import sys
results_seg=sys.argv[1]
annotations_sp1=sys.argv[2]
annotations_sp2=sys.argv[3]
threshold=float(sys.argv[4])
orthologs_file=open('orthologs.txt','w')
partially_orthologs_file=open('partially_orthologs.txt','w')
perfectly_orthologs_file=open('perfectly_orthologs.txt','w')

def ret_dict_geneID_length(annotation_file): #construct a dictionary of the gene id and the length of this 
    annotation_to_split=open(annotation_file,'r')
    asso_GeneID_length={}
    for line in annotation_to_split:
        length=0
        line=line.rstrip()
        field=line.split('\t')
        geneID=field[1]
        exon_start=field[9].split(',')
        exon_end=field[10].split(',')
        exon_start = [x for x in exon_start if x ] #here we exclude the empty line of the list
        exon_end = [x for x in exon_end if x ]
        for i,j in zip(exon_start,exon_end):
            length=length+int(j)-int(i)
        asso_GeneID_length[geneID]=length
    annotation_to_split.close()
    return(asso_GeneID_length)
        
GeneID_length_species1=ret_dict_geneID_length(annotations_sp1)
GeneID_length_species2=ret_dict_geneID_length(annotations_sp2)

def ret_dict_intersections(seg_file): #build a dictionary wich contain the association between the gene ids and the intersections in bp
    results_file=open(seg_file,'r')
    intersections_sp1_2={}
    for line in results_file:
        line=line.rstrip()
        col=line.split('\t')
        compstr=col[7]+'-'+col[9] #concatenation of the 2 corresponding ids
        if compstr not in intersections_sp1_2 :
            intersections_sp1_2[compstr]=int(col[0])
        else: # incrementation with the first field of the seg file wich correspond to the length of the alignment w/o gap
            intersections_sp1_2[compstr]+=int(col[0])
    results_file.close()
    return(intersections_sp1_2)

def calculate_similarity(annot_1,annot_2,intersections,ortho,part_ortho,perf_ortho,thr):
    for id in intersections:
        gene=id.split('-')
        gene1=gene[0]
        gene2=gene[1]
        length_intersections=intersections[id]
        if gene1 in annot_1:
            length_gene1=annot_1[gene1]
            length_gene2=annot_2[gene2]
            unions=length_gene1-length_intersections+length_gene2
            if unions > 0 :
                similarity=float(length_intersections/unions) # here we calculate similarity rate 
                if similarity == 1.0 :
                    perf_ortho.write(gene1+'\t'+gene2+'\t'+str(similarity)+'\n')
                elif similarity < 1 and similarity >= thr :
                    ortho.write(gene1+'\t'+gene2+'\t'+str(similarity)+'\n')
                else :
                    part_ortho.write(gene1+'\t'+gene2+'\t'+str(similarity)+'\n')
            else :
                perf_ortho.write(gene1+'\t'+gene2+'\t1'+'\n')
            #print(f'{gene1} {gene2} {length_intersections} {length_gene1} {length_gene2} {similarity}')
        else :
            length_gene1=annot_2[gene1]
            length_gene2=annot_1[gene2]
            unions=length_gene1-length_intersections+length_gene2
            similarity=float(length_intersections/unions)
            if similarity == 1.0 :
                perf_ortho.write(gene1+'\t'+gene2+'\t'+str(similarity)+'\n')
            elif similarity < 1 and similarity >= thr : # threshold to determine if we have an orthologs or not 
                ortho.write(gene1+'\t'+gene2+'\t'+str(similarity)+'\n')
            else :
                part_ortho.write(gene1+'\t'+gene2+'\t'+str(similarity)+'\n')

 
intersections=ret_dict_intersections(results_seg)

calculate_similarity(GeneID_length_species1,GeneID_length_species2,intersections,orthologs_file,partially_orthologs_file,perfectly_orthologs_file,threshold)
orthologs_file.close()
partially_orthologs_file.close()
perfectly_orthologs_file.close()
#print(unions)
#print(GeneID_length_species1['NM_001011874'])
 