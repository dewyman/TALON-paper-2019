## Run our four GM12878 PacBio datasets through TALON in order to track gene detection by read depth.
Datasets:
- D6
- D7
- D8 (GM12878 Rep 1)
- D9 (GM12878 Rep 2)

Start from database that already contains D8 and D9.
```
cp ../GM12878/full_gencode_v29_2019-03-12.db mult_GM12878.db
```

Then run TALON. 
```
./run_talon.sh
```
