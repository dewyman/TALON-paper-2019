# Follow-up on why the 3' ends of the REST-NRSF reads end earlier than the annotation

## Question 1: Could it be the Classify software?
When we inspected the locaus in question on the genome browser, we noticed that there was a poly-A like sequence in the genome right after the reads ended. This raised the question "could the Classify software have agressively trimmed the reads and removed too much"?

To ask this, I took the HepG2 D4 CCS reads and mapped them to the genome directly. The outcome was that the reads ended in the same place that they do further downstream in the pipeline. So I don't think the bioinformatics are the problem.
