README -- RADIATION 
MLindgren --Oct 2018
---

These are modified radiation scripts written by Steph McAfee back in 2012. Due to some limitations in numpy, we cannot compute above around Latitude 70N. This is due to how their implementation of arccos (acos) is written.  Therefore, I have modified Steph's script and am using it as the DEFAULT way to generate this specific climatology. 

CONSIDER OTHER SCRIPTS (namely Python) in this directory's parent directory as there for legacy only. THIS IS HOW THE RSDS DATA FOR THE CMIP5 RUNS HAVE BEEN CREATED.