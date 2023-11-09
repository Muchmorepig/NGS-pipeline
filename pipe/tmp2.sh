#!/bin/bash

# set the work path
input=${1}
output=${2}
# set the bed file
bed_file=${3}
# set the threads
threads=20

# computerMatirx(This should change)
## scale-regions
computeMatrix scale-regions \
    -S ${input}/*.bw \
    -R ${bed_file} \
    --regionBodyLength 4000 \
    --beforeRegionStartLength 5000 --afterRegionStartLength 5000 \
    --skipZeros --numberOfProcessors ${threads} \
    -o ${output}/matrix_scale_region.scale.gz

## reference-point
computeMatrix reference-point \
    -S ${input}/*.bw \
    -R ${bed_file} \
    --referencePoint TSS \
    -b 5000 -a 5000 \
    --skipZeros --numberOfProcessors ${threads} \
    -o ${output}/matrix_reference_point.reference.gz

# computer multiBigwigSummary
multiBigwigSummary bins -b ${input}/*.bw \
    --numberOfProcessors ${threads} \
    -o ${output}/multibw_results.npz

# plot

## peak plot
plotProfile -m ${output}/matrix_scale_region.scale.gz -out ${output}/scale_region.pdf --perGroup
plotProfile -m ${output}/matrix_scale_region.scale.gz -out ${output}/scale_region_persample.pdf --numPlotsPerRow 4
plotProfile -m ${output}/matrix_reference_point.reference.gz -out ${output}/reference_point_region.pdf --perGroup
plotProfile -m ${output}/matrix_reference_point.reference.gz -out ${output}/reference_point_region_persample.pdf --numPlotsPerRow 4

rm ${output}/matrix_reference_point.reference.gz -rf
rm ${output}/matrix_scale_region.scale.gz -rf

## correlation plot

### plot scatterplot
# plotCorrelation -in ${output}/multibw_results.npz \
#     --corMethod spearman --skipZeros \
#     --whatToPlot scatterplot \
#     --plotTitle "Spearman Correlation" \
#     --removeOutliers \
#     --plotFile ${output}/correlation_spearman_bwscore_scatterplot.pdf

# ### plot heatmap
# plotCorrelation -in ${output}/multibw_results.npz \
#     --corMethod spearman --skipZeros \
#     --whatToPlot heatmap \
#     --plotTitle "Spearman Correlation" \
#     --removeOutliers \
#     --plotNumbers \
#     --plotFile ${output}/correlation_spearman_bwscore_heatmapplot.pdf

rm ${output}/multibw_results.npz -rf
