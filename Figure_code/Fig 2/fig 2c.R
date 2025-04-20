


computeMatrix scale-regions \
              -R gencode.v32.gtf \
              -S OneK1K.sort.bw \
              --metagene \
              -o oneK1K.mat.gz --outFileSortedRegions sortedRegions.bed \
              -b 1000 -a 1000 \
              --sortUsing mean \
              --skipZeros \
              -p 16


plotHeatmap -m oneK1K.mat.gz \
            --samplesLabel Read2 \
            --heatmapWidth 6 --heatmapHeight 9 \
            --colorMap Blues \
            -out oneK1K.pdf
