Scripts for comparing long read and short read splice junctions.

Input:
- High-confidence SJ file from STAR (obtained by mapping short reads to the reference genome)
- SJ file from long reads (achieved by running TranscriptClean utility on TALON-generated GTF)

Output: 
- File listing the number of SJs unique o short reads, unique to long reads, and shared by both. 
