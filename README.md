# PIP-SNP's Pipeline, Source code  and Command Line Run    
Welcome to use **PIP-SNP** to prepprocess your SNP data in the  **Command Line Mode**. In this repository, except this README file, you also should see  a **PDF** file as an user manual which describe the technical detials about **PIP_SNP**   and four sub-directories: **Linux_CB_PIP_SNP_Venue1**, as a **Code::Blocks** project source code file, which can be compiled into a binary  execute command **PIP_SNP_V1** ; **Linux_CB_PIP_SNP_Venue2**, as a **Code::Blocks** project source code file, which can be compiled into a binary   execute command **PIP_SNP_V2** ; **Deep_Sythesis**, as a **Code::Blocks** project source code file, which can be compiled into a binary  execute command **DeepSythesis** ;
**Run_Test**, by which you can play with the three compiled executables,  **PIP_SNP_V1**, **PIP_SNP_V2** and **DeepSythesis** with the example dataset. 
The following are the Linux command lines to run test of the three executables.  
## PIP_SNP_Venue1 
### Step 1.) Get Usage Message
```
$ cd ./Run_Test/PIP_SNP_Venue1_Run_Test
$ ./PIP_SNP_V1 -u
Hello, Dear User! SNP Processing Pipeline for LD Bin Detecting, Missing Genotype Imputing, and or Synthesizing is Starting!
Welcome to use this program to do LD Marker Bin Detecting
The usuage of input parameter arguments are listed as followings:
-u or -U: Output this help usuage message
-g or -G: The full name of Genotype file
-l or -L: The full name of LD Bin Mapping Result file
-i or -I: The Individual number
-r or -R: The Threshold for the pairwise LD R2
-c or -C: the  Correlation Method for a Pairwise Genotypic Markers,0:Pearson_Correlation_R2; 1: LD_D_R2; Default(>2): Pearson_Correlation_R2 
-d or -D: the  Detection Method for a Marker Group(LD Block), 0: Right_Breakthrough; 1: Left_Breakthrough; 2: Left and Right Breakthrough; 3: Left or Right Breakthrough; Default(>3) : Right_Breakthrough
-k or -K: the  K Nearest Neighbor individuals in The KNN method 
-s or -S: the method how to generate the syntesized(binned) genotypic marker, 0: not to synthesize; 1: norm integration of the multiple markers' genotype values; 2: select the optimal one; Default(>2) :norm integration of the multiple markers' genotype values  
-o or -O: the  full name of Imputing results
```
### Step 2.) Run the Command Line
```
 $ ./PIP_SNP_V1 -g ./SNP_Data_Example.txt -l ldmap.txt -o PIP_SNP.txt -i 278 -c 0 -d 0 -r 0.8 -k 10 -s 2
```
## PIP_SNP_Venue2 
### Step 1.) Get Usage Message
```
$ cd ./Run_Test/PIP_SNP_Venue2_Run_Test
$ ./PIP_SNP_V2 -U
Hello, Dear User! Welcome to use SNP Processing Pipeline for LD Bin Filling, Missing Genotype Imputing, and/or Marker Synthesizing!
Welcome to use this program to do LD Marker Bin Detecting
The usuage of input parameter arguments are listed as followings:
-u or -U: Output this help usuage message
-g or -G: The full name of Genotype file
-l or -L: The full name of LD Bin Mapping Result file
-i or -I: The Individual number
-c or -C: the  Correlation Method for a Pairwise Genotypic Markers, 0: Pearson_Correlation_R2; 1: LD_D_R2; Default(>=2): Pearson_Correlation_R2
-k or -K: the  K Nearest Neighbor individuals in The KNN method 
-s or -S: the method how to generate the syntesized(binned) genotypic marker, 0: not to synthesize; 1: norm integration of the multiple markers' genotype values; 2: select the optimal one; default(>2): norm integration   
-o or -O: the  full name of Imputing results
-m or -M: the  full name of LD Mapping results after synthesing
```
### Step 2.) Run the Command Line
```
$ ./PIP_SNP_V2 -g ./SNP_Data_Example.txt -l ./LD_Bin_Map.txt -o PIP_SNP.txt -m ./Synthesis_Map.txt -i 278 -c 0 -k 10 -s 1
```
## DeepSynthesis
###  Step 1.) Get Usage Message
```
$ cd ./Run_Test/Deep_Synthesis_Run_Test
$ ./DeepSynthesis -u
Hello, Dear User! To further  reduce the SNP size, we develop such specific program to deep synthesis SNP Bins by considering the correlationship of a SNP and its neighbor jumped SNP Bins!
Welcome to use this program to do Deep Synthesising
The usuage of input parameter arguments are listed as followings:
-u or -U: Output this help usuage message
-g or -G: The full name of inputing Synthesised SNP file
-l or -L: The full name of inputting  LD Bin file Mapping Synthesising SNPs  
-i or -I: The Individual number
-c or -C: The Method for the pairwise R2 Correlation, 0:Pearson_Correlation_R2; 1: LD_D_R2; Default(>2): Pearson_Correlation_R2 
-r or -R: The Threshold for the pairwise R2 Threshold, a float value (0~1.0)
-s or -S: The Method for the multiple SNP bin Synthesising, 1: norm integration of the multiple markers' genotype values; 2: select the optimal one; Default(0, 0r >2) :norm integration of the multiple markers' genotype values  
-o or -O: The full name of Deep Synthesised SNP results
-m or -M: The full LD Bin mapping file name after deep synthesising
```
 ### Step 2.) Run the Command Line
 ```
$ ./DeepSynthesis -g ./Propressed_SNP.csv -l ./LD_Map.csv -i 132 -c 0 -r 0.2 -s 2 -o ./DS_SNP.txt -m ./DS_Map.txt
 ```
 
If you have any questions, comments, or feedback, please feel free to contact with our developer's team through **pzhaoAT noble DOT org**.  

