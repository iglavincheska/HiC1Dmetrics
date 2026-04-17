python /Users/iglavinchesk/packages/3D-genome-tools/cool2hic.py -i ~/WP2/5_HiC_preevol/2_HiC/data/bin3kb/IPO323_WT_bin3kb_norm.cool -r 3000 -o IPO323_WT_matrix_raw.txt

java -jar /Users/iglavinchesk/packages/3D-genome-tools/juicer_tools.2.20.00.jar pre -r 3000\  /Users/iglavinchesk/packages/3D-genome-tools/IPO323_WT_matrix_KR.txt \  /Users/iglavinchesk/packages/3D-genome-tools/IPO323_WT_bin3kb_norm_KR.hic \
  /Users/iglavinchesk/packages/3D-genome-tools/IPO323_genome.txt

/Users/iglavinchesk/WP2/HiC1Dmetrics/h1d/InteractionFreq.sh /Users/iglavinchesk/packages/3D-genome-tools/juicer_tools.2.20.00.jar /Users/iglavinchesk/packages/3D-genome-tools/IPO323dim2_bin3kb_norm_ICE.hic 4 3000 /Users/iglavinchesk/packages/3D-genome-tools/IPO323_genome.txt /Users/iglavinchesk/WP2/5_HiC_preevol/2_HiC/output/IPO323+dim2_chr4_IF