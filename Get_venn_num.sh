#!/bin/bash
#Fig.1B was made using mergePeaks of HOMER.

ctrl_SE="path to your superEnhancers.txt of anti-IgM 000min data obtained using finPeaks of HOMER"
IgM_SE="path to your superEnhancers.txt of anti-IgM 060min data obtained using finPeaks of HOMER"
ctrl_TE="path to your typicalEnhancers.txt of anti-IgM 000min data obtained using finPeaks of HOMER"
IgM_TE="path to your typicalEnhancers.txt of anti-IgM 060min data obtained using finPeaks of HOMER"

mergePeaks -venn venn_SE.txt -prefix SE ${ctrl_SE} ${IgM_SE}
mergePeaks -venn venn_TE.txt -prefix TE ${ctrl_TE} ${IgM_TE}