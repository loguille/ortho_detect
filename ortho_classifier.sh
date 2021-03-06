# This program require one maf file and two annotations from UCSC genome browser in refgene format.

path=$(dirname $0)
threshold=0.5

while [ -n "$1" ]; do 
    case "$1" in 
    -file) 
        file_org1=$2
        file_org2=$3
        count=0
        for i in $2 $3;
        do 
            count=$(($count+1))
            organism=$(echo $i | cut -d_ -f2 | cut -d. -f1)
            dir=$(dirname $i)
            file=$(basename $i)
            filename=$(basename -- "$i")
            extension="${filename##*.}"
            filename="${filename%.*}"
            if [ $extension == "txt" ]
            then
                echo '############################Importing refgene file to seg format########################################'
                seg-import -c genePred $dir/$file > $dir/$filename.seg
                python3 $path/modify_file.py $dir/$filename.seg $dir/modified_$filename.seg $organism
                seg-sort $dir/modified_$filename.seg > $dir/sort_$filename.seg
                rm $dir/$filename.seg $dir/modified_$filename.seg
            else
                echo '############################Importing GTF file to seg format########################################'
                organism=$(echo $file | cut -d. -f1)
                ./$path/gtfToGenePred $i $dir/tmp.txt
                python3 $path/add_columns.py $dir/tmp.txt $dir/refGene_$organism.txt
                if [ $count -eq 1 ] 
                then
                    file_org1=$dir/refGene_$organism.txt
                else
                    file_org2=$dir/refGene_$organism.txt
                fi
                seg-import -c genePred $dir/refGene_$organism.txt > $dir/refGene_$organism.seg
                python3 $path/modify_file.py $dir/refGene_$organism.seg $dir/modified_refGene_$organism.seg $organism
                seg-sort $dir/modified_refGene_$organism.seg > $dir/sort_refGene_$organism.seg
                rm $dir/refGene_$organism.seg $dir/modified_refGene_$organism.seg $dir/tmp.txt
            echo "done"
            fi
        done;;
    -maf)
        maf_file="$2"
        end_file=$(echo $maf_file| awk -F. '{ print $NF }')
        if [ $end_file == "maf" ];
        then 
            dir_maf=$(dirname $maf_file)
            file=$(basename $maf_file)
            filename=$(basename -- "$maf_file")
            extension="${filename##*.}"
            filename="${filename%.*}"
            echo '############################Importing maf file to seg format########################################'
            seg-import maf $dir_maf/$file > $dir_maf/$filename.seg
            #python3 modify_file_maf.py $dir_maf/$filename.seg $dir_maf/modified_$filename.seg hg38 mm10 #normaly not needed but sometimes maf didn't contain the name of the organism
            #mv $dir_maf/modified_$filename.seg $dir_maf/$filename.seg
            seg-sort $dir_maf/$filename.seg > $dir_maf/sort_$filename.seg
            final_path=$dir_maf/sort_$filename.seg
            rm $dir_maf/$filename.seg
            echo "done"
        else
            echo "it's not the good format"
            exit 0;
        fi;;
    -t) 
        threshold="$2"
    esac
    shift
done;

file_name_1=$(basename $file_org1)
filename_1="${file_name_1%.*}"
file_name_2=$(basename $file_org2)
filename_2="${file_name_2%.*}"
path_1=$(dirname $file_org1)
path_2=$(dirname $file_org2)
path_maf=$(dirname $final_path)
file_maf=$(basename $final_path)
echo "########################### Joining file ###############################"
seg-join $final_path $path_1/sort_$filename_1.seg| seg-swap | seg-sort > intersect.seg
seg-join intersect.seg $path_2/sort_$filename_2.seg | seg-swap | seg-sort > all_intersect.seg
rm $path_1/sort_$filename_1.seg $final_path $path_2/sort_$filename_2.seg intersect.seg
echo "done"
echo "########################### Calculating homology ###############################"
python3 $path/count_intersections_unions.py all_intersect.seg $file_org1 $file_org2 $threshold
echo "done"