ESPDHot is an ensemble machine learning method for protein-DNA binding hotspot prediction.
The authors are Lianci Tao, Tong Zhou, Zhixiang Wu, Fangrui Hu, Shuang Yang, Xiaotian Kong and Chunhua Li.
The performance process includes feature extraction, feature integration and prediction.

In the following, we take a protein-DNA complex (PDB ID: 4CH1) as an example to show the prediction process. For easy description, we give the indications as follows.
'4CH1_A.pdb': the pdb file only containing the protein A chain with the bound DNA removed;
'4CH1_A_bound.pdb': the pdb file containing the information about the protein A chain and its target DNA.

ESPDHot uses the following dependencies:
* python 3.8
* Numpy 1.19.5
* pandas 1.1.5
* joblib
* re
* mlxtend 0.9.1
* xgboost 1.2.1
* sikit-learn 0.24.1
* R 3.6.3
* Matlab R2018b

Step 1: feature extraction

1. Physicochemical characteristics of amino acids
Run ¡®Physicochemical_characteristics.py' (python ./Physicochemical_characteristics.py ./4CH1_A.pdb) with "4CH1_A.pdb" as input file (keep 'Physicochemical_characteristics.py' and '4CH1_A.pdb' in the same folder), to get '4CH1_A_PC.csv'.

2. PSSM
Go to the website 'https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome'. Run PSI-BLAST (select 'PSI-BLAST' for 'Program Selection' item, set E-value = 0.001 and use default parameters for other items) with the sequence of '4CH1_A.pdb' as input to get a PSSM file and rename it as '4CH1_A.asn_matrix.txt'.

3. Coevolutionary
a. Go to the website 'https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome'. Run PSI-BLAST (select 'PSI-BLAST' for 'Program Selection' item with default parameters for other items) with the sequence of '4CH1_A.pdb' as input, then run 'multiple alignment' item with the searched sequences of Query Cover >75%, and E-value <0.0001 as input, and finally download the MSA file and rename it as '4CH1.fa'.
b. run the program 'DCA.m¡¯proposed by Faruck Morcos [(dca(inputfile, outputfile), with inputfile being '4CH1.fa', and outputfile being the output result, 4CH1_A.txt, which should be created in advance].
c. '4CH1_A.txt' needs to be processed by the 'coevolution.py' (python ./coevolution.py ./4CH1_A.txt) to get the '4CH1_A_coevolution.csv' file.

4. Secondary structure
a. Go to the website 'http://sparks-lab.org/server/spider3/'. Run spider3 with the sequence of '4CH1_A.pdb' as input to get the output file and rename it as '4CH1_A_structure.spd33'.
b. Run the python file 'format_transform.py' with '4CH1_A_structure.spd33' as input (python ./format_transform.py ./4CH1_A_structure.spd33) to convert the '4CH1_A_structure.spd33' to '4CH1_A_structure.csv'. 

5. CX/DPX
Go to the website 'https://sourceforge.net/projects/psaia/' to download and install the program 'psaia.exe'.
a. Run 'psaia.exe';
b. Step by step, click 'Structure Analyser', then find 'Analysis Types' and check both 'Analyse as Bound' and 'Analyse by Chain'. All parameters are set to default.
c. Input the pdb file '4CH1_A_bound.pdb' (hydrogen atoms removed) to the program  'psaia.exe', and click 'run' to get the result and rename it as '4CH1_bound.tbl'.
d. Input the pdb file '4CH1_A.pdb' to the program, and click 'run' to get the file and rename it as '4CH1_unbound.tbl'.

6. ASA
a. Go to the website 'http://www.bioinf.manchester.ac.uk/naccess/ ' to download and install the 'Naccess V2.1.1' program on Linux.
b. Run 'Naccess V2.1.1' on Linux and input corresponding pdb file (e.g. 4CH1_A_bound.pdb) with the command:'./naccess 4CH1_A_bound.pdb' to get the .rsa format result file and renamed it as '4ch1_A_bound.rsa'.
c. Repeat step b and get another file '4ch1_A_unbound.rsa' through inputting the '4CH1_A.pdb' file.

7. Solvent exposure
Run the website 'https://sunflower.kuicr.kyoto-u.ac.jp/~sjn/hse/webserver.html' by inputting the sequence of '4CH1_A.pdb' (select 'PSI-BLAST+PSIPRED+AA+W+L' for 'SVR Models' item) to get the result from the e-mail that you have submitted and name the result file as '4CH1_A_HSE.txt'.

8. Interface preference
a. Run the website 'https://swift.cmbi.umcn.nl/gv/dssp/index.html' to get the dssp file of the 4CH1 and only reserve the structure information of the A chain of the protein and rename it as '4CH1_A.dssp'
b. Run IP.py (python ./IP.py ./4CH1_A.dssp ./60_4P.data) with 4CH1_A.dssp and 60_4P.data as input file (keep 'IP.py', '4CH1_A.dssp' and '60_4P.data' in the same folder) to get¡®4CH1_A_IP.csv¡¯.

9. Network topological feature
Go to the website 'http://www.sysbio.org.cn/NACEN/' to download 'dssp-3.0.0.exe' and R package 'NACEN', and see the instruction on the website to perform the following operations.
a. Run R, install the package 'NACEN' and then library(NACEN);
b. Use function 'NACENConstructor' to construct NACEN with the parameter "WeightType = 'SAS'" (node-weighted AAN based on residue solvent accessibility);
c. Run function 'NACEAnanlyzer' to obtain a result file named it as '4CH1_A_SAS.txt';
d. Repeat the above procedures with "WeightType = 'SAS'" replaced by "WeightType = 'Hydrophobicity'", "WeightType = 'Polarity'" and "WeightType = 'Mass'" to obtain the result files named they as '4CH1_A_Hydrophobicity.txt', '4CH1_A_Polarity.txt' and '4CH1_A_Mass.txt', respectively.

10. Dynamics
Run pfGNM.py (python ./pfGNM.py ./4CH1_A.pdb) to get a result '4CH1_A_dynamics.csv'

Step 2: prediction integration

After all features are extracted, all the feature files obtained above including '4CH1_A_PC.csv', '4CH1_A.asn_matrix.txt', '4CH1_A_coevolution.csv', '4CH1_A_structure.csv', '4CH1_bound.tbl', '4CH1_unbound.tbl', '4ch1_A_bound.rsa', '4ch1_A_unbound.rsa', '4CH1_A_HSE.txt', '4CH1_A_IP.csv', '4CH1_A_SAS.txt', '4CH1_A_Hydrophobicity.txt', '4CH1_A_Polarity.txt', '4CH1_A_Mass.txt', '4CH1_A_dynamics.csv'.
Put all feature files mentioned above in the 'Input' folder where 'Feature_integration.py' is located.
Run Feature_integration.py (python ./Feature_integration.py ./4CH1_A.pdb) with the 'Input' folder and '4CH1_A.pdb' in the same folder where 'Feature_integration.py' is located to get the '4CH1_A_feature.csv' file in the Output file. 

Step 3: prediction

Put '4CH1_A_feature.csv' in the 'Input_data' folder, and keep 'predict.py', 'Model_file' folder and 'Input_data' folder in the same path.
Run the Command: 'python ./predict.py ./Input_data/4CH1_A_feature.csv' to get the output "./Result/4CH1_A_hotspot_result.csv" where the predicted hotspots are stored.

It should be pointed out that although ESPDHot can give the prediction results for all the residues, only the results corresponding to the interfacial residues are significant. 

Help
For any questions, please contact us at chunhuali@bjut.edu.cn.
