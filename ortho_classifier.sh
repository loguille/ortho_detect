if [ $# -ne 3 ]
then
    echo "Invalid number of argument (must be 3)"
else
    sh compile_ortho_decider.sh $1 $2
    organism=$(echo $3 | cut -d_ -f2 | cut -d. -f1)
    dir=$(dirname $3)
    filename_ref=$(basename -- "$3")
    extension="${filename_ref##*.}"
    filename_ref="${filename_ref%.*}"
    file_seg="$filename_ref.seg"
    echo "#############################Importing RefGene file to seg format#####################################"
    seg-import genePred $dir/$filename_ref.$extension > $dir/$file_seg
    newfilename='modified_'$file_seg
    python3 modify_file.py $dir/$file_seg $dir/$newfilename $organism
    seg-sort $dir/$newfilename > $dir/sort_$file_seg
    rm $dir/$newfilename $dir/$file_seg
    echo "########################### Joining file ###############################"
    seg-join sort_swap_intersections_maf_refgene.seg $dir/sort_$file_seg > intersect.seg
    seg-swap intersect.seg > tmp.seg
    seg-sort tmp.seg > all_intersect.seg
    rm sort_swap_intersections_maf_refgene.seg $dir/sort_$file_seg intersect.seg tmp.seg
    echo "###################################### Decide orthology ###################################"
    python3 position_comparison.py all_intersect.seg $2 $3 
fi