# Performance testing the TALON pipeline

1. See input_files directory README for instructions to generate the input files.

2. Generate and run Minimap2 scripts for the variously size inputs
```
FPATH=/pub/dwyman/TALON-paper-2019/performance_testing

python ${FPATH}/generate_Minimap2_job.py --dir ${FPATH}/input_files --input 1000 --logs ${FPATH}/logs --outdir ${FPATH}/runs
```
