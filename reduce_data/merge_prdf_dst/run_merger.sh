#!/bin/bash

#remove the echo
#merge the dst tree dump with the full PRDF dump file, edit as neccessary
./merge.exe ./359711_DST_dump.txt ../../data/run_12/prdf_data/359711_merged_corrected.txt ./merged_trees/359711_dst_prdf.txt
./merge.exe ./360879_DST_dump.txt ../../data/run_12/prdf_data/360879_merged_corrected.txt ./merged_trees/360879_dst_prdf.txt
./merge.exe ./362492_DST_dump.txt ../../data/run_12/prdf_data/362492_merged_corrected.txt ./merged_trees/362492_dst_prdf.txt
./merge.exe ./364636_DST_dump.txt ../../data/run_12/prdf_data/364636_merged_corrected.txt ./merged_trees/364636_dst_prdf.txt
./merge.exe ./365866_DST_dump.txt ../../data/run_12/prdf_data/365866_merged_corrected.txt ./merged_trees/365866_dst_prdf.txt
./merge.exe ./366605_DST_dump.txt ../../data/run_12/prdf_data/366605_merged_corrected.txt ./merged_trees/366605_dst_prdf.txt
./merge.exe ./367138_DST_dump.txt ../../data/run_12/prdf_data/367138_merged_corrected.txt ./merged_trees/367138_dst_prdf.txt
