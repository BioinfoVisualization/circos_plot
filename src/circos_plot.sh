#!/bin/sh
# ============================================================================ #
# circos_plot.sh                                                               #
# Circos plot wrapper. Calls the python script adj4circos.py to generate a     #
# usage adjacency matrix between V and J usages. It then calls the R script    #
# circos_plot.R to generate the plot and export it as a pdf figure.            #
# Parameters:                                                                  #
#     chain : string. The letter of the chain (A,B,G or D)                     #
#     patient : int. Patient number (0-7). If zero, no patient is              #
#               selected.                                                      #
#     tissue : string. The name of the tissue (MUSL, PB, both)
#     TC : int. Threshold for sample filtering. If zero no filtering is done.  #
#                                                                              #
# Author: Juan Sebastian Diaz Boada                                            #
# Creation Date: 4/11/2021                                                     #
# ============================================================================ #
chain=$1
patient=$2
tissue=$3
TC=$4
if [ $patient -eq 0 ]
then
  echo "No patient specified. Using complete dataset..."
  pat=""
else
  echo "Patient identified. Subsetting dataset..."
  pat="_sc${patient}"
fi
if [ $tissue = "both" ]
then
  echo "No tissue specified. Using complete dataset..."
  tis=""
else
  echo "Tissue identified. Subsetting dataset..."
  tis="_${tissue}"
fi
TCR="../data/output_data/TCR_metadata.tsv"
if [ $TC -eq 0 ]
then
  echo "No threshold. Not performing any sample filtering..."
  ./adj4circos.py  $TCR $chain --patient $patient --tissue $tissue --no-filter
  in_file="../data/output_data/circos_adj/circos_adjacency_${chain}${pat}${tis}.csv"
  out_file="../figures/circos/circos_${chain}${pat}${tis}.pdf"
else
  echo "Threshold identified. Filtering..."
  ./adj4circos.py  $TCR $chain --patient $patient --tissue $tissue --filter --thresh $TC
  in_file="../data/output_data/circos_adj/circos_adjacency_${chain}${pat}${tis}_thresh_${TC}.csv"
  out_file="../figures/circos/circos_${chain}${pat}${tis}_thresh_${TC}.pdf"
fi
Rscript --vanilla circos_plot.R $in_file $out_file
