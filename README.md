# OffsetTreeMaker
L1 jet corrections derived using root tuples
Step 1: Making tuples with OffsetTreeMaker.cc
     for Data:
     a) extract corresponding json pileup from this site:
     https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/
     b) Convert json pileup format to ASCII using JSONtoASCII.py "pileup_latest.txt"
     c) modify plugins/parsePileUpJSON2.h with corresponding pileup json accordingly.
     d) add global tags to run_offset.py
     e) process entire sample using crab_run_offset.py
     for MC:
     repeat above with skipping steps a), b), c)
     
     
step 2: after making ntuples we process them using histomaker.cc:
    for data:
    nohup histomaker false 0.4(0.8) data.root
    for MC:
    nohup histomaker true 0.4(0.8) data.root mc.root (data.root is the file that you reweight wrt)
    
step3:
    pick a range of pileup to process samples:
    e.g. (10,30)
    
    a) root -l -b -q 'offsetpT.c (10,30,30,30)'
      (pick all and chs)
    b) root -l -b -q 'scalefactor.c ("all",10,30)' and root -l -b -q 'scalefactor.c ("chs",10,30)' 
    
    c) root -l -b -q 'l1fastjet.c("all")' and root -l -b -q 'l1fastjet.c("chs")'
