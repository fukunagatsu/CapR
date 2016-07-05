# CapR
CapR calculates probabilities that each RNA base position is located within each secondary structural context for long RNA sequences. Six categories of RNA structures were taken into account, stem part (S), hairpin loop (H), Bulge loop (B), Internal loop (I), multibranch loop (M), and exterior loop (E). We defined a structural profile of an RNA base by a set of six probabilities that the base belongs to each category. 

##Version
Version 1.1.1 (2016/07/05)

##Acknowledgements
We used a portion of the source code from the Vienna RNA package 1.8.5. We thank Dr. Ivo L. Hofacker, Vienna RNA package development group, and the Institute for Theoretical Chemistry of the University of Vienna. You can download the source code for Vienna RNA Package from http://www.tbi.univie.ac.at/RNA/

##Usage
    ./CapR <input_file> <output_file> <maximal_span>

input_file : input RNA sequences in fasta format.  
output_file : output file name.  
maximal_span  : maximal length between bases that form base pairs.  

CapR output is a space delimited format with a fasta-like header line.  
The probabilities (prob1 prob2 ...) are order by the sequence positions.  

Example:  
\>M33000.1/55-110  
Bulge 0 3.6242e-10 2.92873e-07 7.63119e-06 8.06009e-05 ...  
Exterior 1 0.108771 0.000237711 0.00018179 0.00010304 ...  
Hairpin 0 2.122e-08 8.20875e-08 1.69936e-07 3.61726e-07 ...  
Multibranch 0 4.48371e-12 3.45769e-05 8.38716e-05 0.0165036 ...  
Stem 2.23919e-07 0.891229 0.999688 0.999701 0.98264 ...  

##License
This software is released under the MIT License, see LICENSE.txt.

##Changelogs
2016/07/05 Version 1.1.1 bug fix: I fixed a bug in "fastafile_reader.cpp"  
2016/06/30 Version 1.1.0 was released.


## Reference
Tsukasa Fukunaga, Haruka Ozaki, Goro Terai, Kiyoshi Asai, Wataru Iwasaki, and Hisanori Kiryu. "CapR: revealing structural specificities of RNA-binding protein target recognition using CLIP-seq data." *Genome Biology*, **15**, R16. (2014)
