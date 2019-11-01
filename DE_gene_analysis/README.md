## Differential Gene identification between different sequencing platforms
```bash
mkfir figures
```

0. Reformat tables (once you get them from Gaby)
```bash
pb=pb_GM12878_K562.tsv
ont=ont_GM12878_K562.tsv
ill=illumina_GM12878_K562.tsv

cp et_PacBio_GM12878vsK562.txt $pb
cp et_ONT_GM12878vsK562 $ont 
cp et_illumina_GM12878vsK562.txt $ill

# add 'gene' column
printf "gene " > temp 
cat $pb >> temp 
mv temp $pb
printf "gene " > temp 
cat $ont >> temp 
mv temp $ont
printf "gene " > temp 
cat $ill >> temp 
mv temp $ill

# remove all quotations
sed -i temp 's/"//g' $pb
sed -i temp 's/"//g' $ont
sed -i temp 's/"//g' $ill

# replace all tabs or spaces with ,
sed -i temp 's/ /,/g' $pb
sed -i temp 's/ /,/g' $ont
sed -i temp 's/ /,/g' $ill
```

1. Plot differentially expressed genes between HepG2 and K562
```bash
python ma_plot_DGs.py \
	-pb pb_HepG2_K562_DEG.tsv \
	-ont ont_HepG2_K562_DEG.tsv \
	-ill illumina_HepG2_K562_DEG.tsv \
	-sample "HepG2 to K562"
```

Plot differentially expressed genes between GM12878 and K562
```bash
python ma_plot_DGs.py \
	-pb $pb \
	-ont $ont \
	-ill $ill \
	-sample "GM12878 to K562"
```

```python
	# MA plot colored by which technologies identified each gene as DE
	plt.figure(figsize=(8.5,8.5))
	ax = sns.scatterplot(data=df, x='log2FC_a', y='log2CPM_a', hue='de_category',
		 	alpha=0.5, size=3, edgecolor=None)
	plt.xlabel('{} log2-fold change in {}'.format(sample_name, namea))
	plt.ylabel('log(Counts per million) in {}'.format(namea))
	handles, labels = ax.get_legend_handles_labels()
	ax.legend(handles=handles[1:], labels=labels[1:])
	# fig = ax.get_figure()
	plt.savefig('test_pb.png')
	exit()
```
 <!-- <img align="center" width="250" src="figures/GM12878_venn.png"> <img align="center" width="250" src="figures/HepG2_venn.png">  -->