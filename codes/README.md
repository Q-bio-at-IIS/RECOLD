REpertoire COmparison in Low Dimensions (RECOLD)
===================


RECOLD is a method to compare T cell receptor (TCR) repertoires among different samples.   
Here, we release the accompanying python code for the paper: R. Yokota, Y. Kaminaga, and T.J. Kobayashi, "Quantification of inter-sample differences in T cell receptor repertoires using sequence-based information", *Front. Immunol.*, 8:1500; doi: 10.3389/fimmu.2017.01500.

----------

### Installation instructions

Our code requires the sequence alignment library "parasail" (https://github.com/jeffdaily/parasail).  
Before running the code, please install the library by following the instructions on the parasail Github page on the above link. 
Our software consists with the following directory tree;

```
.
└── codes
    ├── Data
    │   └── README.txt
    ├── MatData
    ├── main.py
    └── mymodule
        ├── mymodules.py
        └── mymodules.pyc
```
In our article,  we used the published data from the other papers.   To confirm that our code works well, please download the supplementary data1 of Rempala et al., (2011)[^1] from the website of Journal of Theoretical Biology, and storage the downloaded csv file in the "Data" directory.
When you run the main.py,  the output files are generated in the MatData directory.  

### Please cite RECOLD as:
Yokota R, Kaminaga Y and Kobayashi TJ (2017) Quantification of Inter-Sample Differences in T-Cell Receptor Repertoires Using Sequence-Based Information. Front. Immunol. 8:1500. doi: 10.3389/fimmu.2017.01500

### Reference
 [^1]: Rempala GA, Seweryn M, Ignatowicz ., “Model for comparative analysis of antigen receptor repertoires.”, Journal of Theoretical Biology, Vol.269., pp.1-15, (2011).




