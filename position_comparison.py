import sys
#python3 position_comparison.py head_all_intersect_hg38mm10_simplified.seg  HumanvsMouse/starting-file/refGene_Human.txt HumanvsMouse/starting-file/refGene_Mouse.txt

orthologs_file=open('orthologs.txt','w')
partially_orthologs_file=open('partially_orthologs.txt','w')

def make_dictionary_id_exons_junctions(file_annot):
    file_to_split=open(file_annot,'r')
    dictionary_refgene={}
    for line in file_to_split:
        elem=line.split('\t')
        position_start=elem[9].split(',')
        position_end=elem[10].split(',')
        position_start = [x for x in position_start if x ]
        position_end = [x for x in position_end if x ]
        position=(elem[2],elem[3],position_start,position_end)
        if elem[1] not in dictionary_refgene:
            dictionary_refgene[elem[1]]=position
        #else :
        #    position=list(position)+list(dictionary_refgene[elem[1]])  bien que certains id soient en plusieurs nombres ils correspondent a des transcrits different
        #    dictionary_refgene[elem[1]]=position
    file_to_split.close()   
    return(dictionary_refgene)

def make_list_with_seg_file(file_seg):
    seg_file=open(file_seg,'r')
    list_seg_file=[]
    for line in seg_file:
        line=line.rstrip() 
        elem=line.split('\t')
        list_seg_file.append(elem)
    seg_file.close()
    return(list_seg_file)


species1_annotations=make_dictionary_id_exons_junctions(sys.argv[2])#contient les informations sur les annotations 
species2_annotations=make_dictionary_id_exons_junctions(sys.argv[3])
seg_list=make_list_with_seg_file(sys.argv[1]) 

def make_dictionnary_from_seg_list(list_seg):
    simple_dictionnary={}
    for elem in list_seg: 
        if elem[9] not in simple_dictionnary:
            simple_dictionnary[elem[9]]=[elem]
        else:
            simple_dictionnary[elem[9]].append(elem)
    return(simple_dictionnary)

def construct_dictionnary_id_id(list_seg):
    dico={}
    for elem in list_seg:
        if elem[7] not in dico:
            dico[elem[7]]=[elem[9]]
        else:
            dico[elem[7]].append(elem[9])
    for key in dico:
        dico[key]=list(set(dico[key]))
    return(dico)


def make_dictionary_with_seg_file(list_seg): 
    dictionary_segfile={}
    for elem in list_seg:
        if elem[7] not in dictionary_segfile:
            dictionary_segfile[elem[7]]=1
        else:
            dictionary_segfile[elem[7]]+=1
    return(dictionary_segfile)



id_nb_line=make_dictionary_with_seg_file(seg_list) #create a dictionary containing the id and the number of line associated to this id 
association_id_line=make_dictionnary_from_seg_list(seg_list)#create a dictionnary where each id of the 2nd species contain all the lines associated
association_id_id=construct_dictionnary_id_id(seg_list)# create a dictionnary wher each gene id of the 1st species is associated to all the gene id of the 2nd species which it have lines with 

#print(association_id_id)

#elem[0]    elem[1] elem[2]     elem[3]         elem[4]         elem[5] elem[6]     elem[7]     elem[8]     elem[9]     elem[10]
#897     hg38.chr1       69111   mm10.chr2       111483431       8       583     NM_001005484    21      NM_146404       0 (seg_file_with_annotations)
def make_decision(dict_id_id, dict_id_line,dict_annot_sp1,dict_annot_sp2):
    #cette partie expose comment acceder aux differents element 
    for id in dict_id_id:
        print('############################# new information #############################')
        values=dict_id_id[id]
        print(id)
        print(dict_annot_sp1[id]) #('chr1', '+', ['69090'], ['70008'])
        for i in values:
            print(dict_id_line[i])#[['897', 'hg38.chr1', '69111', 'mm10.chr2', '111483431', '8', '583', 'NM_001005484', '21', 'NM_146404', '0']]
            print(dict_annot_sp2[i])#('chr2', '+', ['111483431'], ['111484328'])

#make_decision(association_id_id,association_id_line,species1_annotations,species2_annotations)

def make_decision1vs1(dict_id_id, dict_id_line,dict_annot_sp1,dict_annot_sp2,id_gene1,id_gene2):
    #this method works for 1 to 1 exons
    results='non-orthologous'
    start_gene_1=int(dict_annot_sp1[id_gene1][2][0])
    end_gene_1=int(dict_annot_sp1[id_gene1][3][0])
    start_gene_2=int(dict_annot_sp2[id_gene2][2][0])
    end_gene_2=int(dict_annot_sp2[id_gene2][3][0])
    for line in dict_id_line[id_gene2]:
        start_cmp_species1=abs(int(line[2]))
        end_cmp_species1=abs(int(line[2])+int(line[0]))
        start_cmp_species2=abs(int(line[4]))
        end_cmp_species2=abs(int(line[4])+int(line[0]))
        if start_cmp_species1 <= start_gene_1 and start_cmp_species2 <= start_gene_2 and end_cmp_species1 >= end_gene_1 and end_cmp_species2 >= end_gene_2:
            results='orthologous'
            break     
        #elif start_cmp_species1 > end_gene_1 or end_cmp_species1 < start_gene_1 or start_cmp_species2 > end_gene_2 or end_cmp_species2 <start_gene_2 :
        #    results='non-orthologous'      
        else :
            results='partially_orthologous'
    return(results)          
"""        
    if results == 'orthologous':
        return(results)
        #print(f'Results between {id_gene1} and {id_gene2} \n This genes are orthologous ################################################')
    else :
        return(results)
        #print(f'Results between {id_gene1} and {id_gene2} \n This genes are non orthologous')
"""



def make_decision_multvsmult(dict_id_id, dict_id_line,dict_annot_sp1,dict_annot_sp2,id_gene1,id_gene2): 
    liste_start_exons1=list(dict_annot_sp1[id_gene1][2])
    liste_end_exons1=list(dict_annot_sp1[id_gene1][3])
    liste_start_exons2=list(dict_annot_sp2[id_gene2][2])
    liste_end_exons2=list(dict_annot_sp2[id_gene2][3])
    for line in dict_id_line[id_gene2]:
        start_cmp_species1=abs(int(line[2]))
        end_cmp_species1=abs(int(line[2])+int(line[0]))
        start_cmp_species2=abs(int(line[4]))
        end_cmp_species2=abs(int(line[4])+int(line[0]))
        for start_exons1,end_exons1 in zip(liste_start_exons1,liste_end_exons1):
            for start_exons2, end_exons2 in zip(liste_start_exons2,liste_end_exons2):
                #print(f'debut exons gene1 {start_exons1} fin exons gene1 {end_exons1} debut comp especes 1 {start_cmp_species1} fin cmp especes 1 {end_cmp_species1} deb comp esp 2 {start_cmp_species2} fin cmp esp2 {end_cmp_species2} debut exons gene2 {start_exons2} fin exons gene2 {end_exons2}')
                if start_cmp_species1 <= int(start_exons1) and end_cmp_species1 >= int(end_exons1) and start_cmp_species2 <= int(start_exons2) and end_cmp_species2 >= int(end_exons2):
                    #print(f'{liste_start_exons1} match : {start_exons1}')
                    liste_start_exons1.remove(start_exons1)
                    liste_end_exons1.remove(end_exons1)
                    liste_start_exons2.remove(start_exons2)
                    liste_end_exons2.remove(end_exons2)
                    #print(f' nouvelle liste : {liste_start_exons1}')
    if len(liste_start_exons1) == 0 or len(liste_start_exons2) == 0:
        return('orthologous')
        #print(f'Results between {id_gene1} and {id_gene2} \n This genes are orthologous {len(liste_start_exons1)} {len(liste_start_exons2)} numbers of exons {len(dict_annot_sp1[id_gene1][2])} {len(dict_annot_sp2[id_gene2][2])} #####################')
    else:
        return('partially_orthologous')
        #print(f'Results between {id_gene1} and {id_gene2} \n This genes are partially orthologous\n {len(liste_start_exons1)} {len(liste_start_exons2)} numbers of exons {len(dict_annot_sp1[id_gene1][2])} {len(dict_annot_sp2[id_gene2][2])}') 

        
def make_decision_same_mumber(dict_id_id, dict_id_line,dict_annot_sp1,dict_annot_sp2,id_gene1,id_gene2):
    liste_start_exons1=list(dict_annot_sp1[id_gene1][2])
    liste_end_exons1=list(dict_annot_sp1[id_gene1][3])
    liste_start_exons2=list(dict_annot_sp2[id_gene2][2])
    liste_end_exons2=list(dict_annot_sp2[id_gene2][3])
    for line in dict_id_line[id_gene2]:
        start_cmp_species1=abs(int(line[2]))
        end_cmp_species1=abs(int(line[2])+int(line[0]))
        start_cmp_species2=abs(int(line[4]))
        end_cmp_species2=abs(int(line[4])+int(line[0]))
        for start_exon1,end_exon1,start_exon2,end_exon2 in zip(liste_start_exons1,liste_end_exons1,liste_start_exons2,liste_end_exons2):
            #print(f'debut exons gene1 {start_exon1} debut comp especes 1 {start_cmp_species1}\tfin exons gene1 {end_exon1}  fin cmp especes 1 {end_cmp_species1}\tdeb comp esp 2 {start_cmp_species2} debut exons gene2 {start_exon2}\tfin cmp esp2 {end_cmp_species2} fin exons gene2 {end_exon2}')
            if start_cmp_species1 <= int(start_exon1) and end_cmp_species1 >= int(end_exon1) and start_cmp_species2 <= int(start_exon2) and end_cmp_species2 >= int(end_exon2):
                #print(f'debut exons gene1 {start_exon1} debut comp especes 1 {start_cmp_species1}\tfin exons gene1 {end_exon1}  fin cmp especes 1 {end_cmp_species1}\tdeb comp esp 2 {start_cmp_species2} debut exons gene2 {start_exon2}\tfin cmp esp2 {end_cmp_species2} fin exons gene2 {end_exon2}')
                liste_start_exons1.remove(start_exon1)
                liste_end_exons1.remove(end_exon1)
                liste_start_exons2.remove(start_exon2)
                liste_end_exons2.remove(end_exon2)
    if len(liste_start_exons1) == 0 :
        return('orthologous')
        #print(f'Results between {id_gene1} and {id_gene2} \n This genes are orthologous ##################################')
    else :
        return('partially_orthologous')
        #print(f'Results between {id_gene1} and {id_gene2} \n This genes are partially orthologous\n {len(liste_start_exons1)} {len(dict_annot_sp1[id_gene1][2])}')   

for gene in association_id_id:
    associated_gene_list=association_id_id[gene]
    #print(associated_gene_list)
    for associated_gene in associated_gene_list:
        #make_decision_multvsmult(association_id_id,association_id_line,species1_annotations,species2_annotations,gene,associated_gene)
        #print(associated_gene)
        if len(species1_annotations[gene][2]) == 1 and len(species2_annotations[associated_gene][2]) == 1 :
            decision=make_decision1vs1(association_id_id,association_id_line,species1_annotations,species2_annotations,gene,associated_gene)
            if decision=='orthologous':
                orthologs_file.write(gene+'\t'+associated_gene+'\n')
            else:
                partially_orthologs_file.write(gene+'\t'+associated_gene+'\n')
        elif len(species1_annotations[gene][2])!=len(species2_annotations[associated_gene][2]) and (len(species1_annotations[gene][2]) != 1 or len(species2_annotations[associated_gene][2]) != 1):
            decision=make_decision_multvsmult(association_id_id,association_id_line,species1_annotations,species2_annotations,gene,associated_gene)
            if decision=='orthologous':
                orthologs_file.write(gene+'\t'+associated_gene+'\n')
            else:
                partially_orthologs_file.write(gene+'\t'+associated_gene+'\n')
        else : 
            decision=make_decision_same_mumber(association_id_id,association_id_line,species1_annotations,species2_annotations,gene,associated_gene)
            if decision=='orthologous':
                orthologs_file.write(gene+'\t'+associated_gene+'\n')
            else:
                partially_orthologs_file.write(gene+'\t'+associated_gene+'\n')
orthologs_file.close()
partially_orthologs_file.close()