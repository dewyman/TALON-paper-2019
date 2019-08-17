Since we ran TranscriptClean on different SMRT cells separately, it is necessary to concatenate the outputs in order to generate one report.

```
cat ../../../../../pipeline/D8/TC_v1.0.7/3_C01/3_C01_clean.TE.log \
    <(tail -n+2 ../../../../../pipeline/D8/TC_v1.0.7/4_D01/4_D01_clean.TE.log) \
    <(tail -n+2 ../../../../../pipeline/D9/TC_v1.0.7/1_A01/1_A01_clean.TE.log) \
    <(tail -n+2 ../../../../../pipeline/D9/TC_v1.0.7/2_B01/2_B01_clean.TE.log) \
    > GM12878_clean.TE.log

cat ../../../../../pipeline/D8/TC_v1.0.7/3_C01/3_C01_clean.log \
    <(tail -n+2 ../../../../../pipeline/D8/TC_v1.0.7/4_D01/4_D01_clean.log) \
    <(tail -n+2 ../../../../../pipeline/D9/TC_v1.0.7/1_A01/1_A01_clean.log) \
    <(tail -n+2 ../../../../../pipeline/D9/TC_v1.0.7/2_B01/2_B01_clean.log) \
    > GM12878_clean.log

```

Run report script
```
Rscript ~/TranscriptClean-1.0.7/generate_report.R GM12878
```
