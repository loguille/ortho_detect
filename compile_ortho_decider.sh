#### This program take as an entry 2 files the first is a maf file from the 2 organism that has been aligned and the others RefGene format for the organism
## this program use seg suite develop by Pr Martin Frith



if [ "$#" -ne 2 ] ;
then

    echo 'invalid number of arguments (must be 2)';
else
    for file in $1 $2 ; 
    do
        end_file=$(echo $file| awk -F. '{ print $NF }')
        if [ $end_file == "maf" ]; 
        then
            dir_maf=$(dirname $file)
            maf_file=$(basename $file)
            filename=$(basename -- "$file")
            extension="${filename##*.}"
            filename="${filename%.*}"
            echo '############################Importing maf file to seg format########################################'
            seg-import maf $dir_maf/$maf_file > $dir_maf/$filename.seg
        else 
            organism=$(echo $file | cut -d_ -f2 | cut -d. -f1)
            dir=$(dirname $file)
            filename_ref=$(basename -- "$file")
            extension="${filename_ref##*.}"
            filename_ref="${filename_ref%.*}"
            file_seg="$filename_ref.seg"
            echo "#############################Importing RefGene file to seg format#####################################"
            seg-import genePred $dir/$filename_ref.$extension > $dir/$file_seg
            newfilename='modified_'$file_seg
            python3 modify_file.py $dir/$file_seg $dir/$newfilename $organism
            seg-sort $dir/$newfilename > $dir/sort_$file_seg
            rm $dir/$newfilename $dir/$file_seg

        fi;
    done;
    echo "#######################################Joining files####################################################"
    seg-join $dir_maf/$filename.seg $dir/sort_$file_seg > intersections_maf_refgene.seg
    seg-swap intersections_maf_refgene.seg > swap_intersections_maf_refgene.seg
    seg-sort swap_intersections_maf_refgene.seg > sort_swap_intersections_maf_refgene.seg
    rm intersections_maf_refgene.seg swap_intersections_maf_refgene.seg $dir/$filename.seg $dir/sort_$file_seg
     
fi

