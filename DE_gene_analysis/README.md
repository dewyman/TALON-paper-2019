## Differential Gene identification between different sequencing platforms
```bash
mkfir figures
```

1. Plot differentially expressed genes between HepG2 and K562
```bash
python ma_plot_DGs.py \
	-pb pb_HepG2_K562_DEG.tsv \
	-ont ont_HepG2_K562_DEG.tsv \
	-ill illumina_HepG2_K562_DEG.tsv \
	-sample "HepG2 vs. K562"
```

 <!-- <img align="center" width="250" src="figures/GM12878_venn.png"> <img align="center" width="250" src="figures/HepG2_venn.png">  -->