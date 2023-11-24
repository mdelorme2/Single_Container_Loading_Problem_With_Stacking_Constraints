This repository contains the code for all algorithms discussed in the paper "Exact decomposition approaches for a single container loading problem with stacking constraints and medium-sized weakly heterogeneous items" by Maxence Delorme and Joris Wagenaar. 

Our algorithms are coded in C++ and use the commercial solver Gurobi for the ILP models and Cplex for the CP models. 
The code is divided over 17 folders, each containing the code of one method. The different folders correspond to the following methods in our paper:
- 0_NOTHING       		| COMPACT-NOP
- 1_COMPACT_BPP         | COMPACT
- 1_COMPACT_CSP 		| COMPACT-CSP
- 2_SEQ					| SEQ
- 2_SEQ_KEP   		    | SEQ-KEP
- 3_DEC1     			| DEC1-CB
- 3_DEC1_SEQ     		| DEC1
- 4_DEC2				| DEC2-CB
- 4_DEC2_SEQ			| DEC2
- 5_DEC3				| DEC3-CB
- 5_DEC3_SEQ			| DEC3
- 5_DEC3_SEQ_00			| DEC3-0-0
- 5_DEC3_SEQ_01			| DEC3-1-0 
- 5_DEC3_SEQ_10			| DEC3-0-1
- 5_DEC3_SEQ_11			| DEC3-1-1

Each folder contains the same substructure. For example, 2_SEQ contains the following files:
- Allocation.cpp: contains a number of secondary functions (this file is usually the same for each of the 17 main folders)
- Allocation.h: the header file corresponding to Allocation.cpp (this file is usually the same for each of the 17 main folders)
- main.cpp: the front-end code for using the method SEQ
- main.h: the header file corresponding to SEQ.cpp
- time.cpp: contains the functions to compute the running time of the algorithms (this file is the same for each of the 17 main folders)
- time.h: the header file corresponding to time.cpp (this file is the same for each of the 17 main folders)
- makefile: used for compiling under linux (it needs to be updated by the user)

Once compiled, the following command can be used to run the algorithm:
	./PROGRAM "./PATH_INSTANCE" "NAME_INSTANCE" "./PATH_AND_NAME_OUTPUT_GENERAL" "./PATH_AND_NAME_OUTPUT_PICTURE"  
where
- PROGRAM is the name of the compiled software (e.g., SEQ)
- ./PATH_INSTANCE is the relative path of the folder where the instance to solve is located
- NAME_INSTANCE is the name of the instance to solve
- ./PATH_AND_NAME_OUTPUT_GENERAL is the name of the file (together with its relative path) where performance metrics (such as the optimality status, the CPU time required, or the number of variables) are stored after solving an instance
- ./PATH_AND_NAME_OUTPUT_PICTURE is the name of the file (together with its relative path) where the file that serves as input for the picture generator is stored

Moreover, "_INPUT.rar" contains a txt-file for each of our test instances. There are 6 main folders, each corresponding to a different instance type:
- GENERATED_DIM	    (Varying the truck dimensions W x L x H)
- GENERATED_LAY	    (Varying the item packing dimensions w_i and l_i range)
- GENERATED_MUL	    (Varying the number of item types m)
- GENERATED_QUA	    (Varying the number of items n)
- REALISTIC 	    (Real instances)
- REALISTIC_C	    (Real instances with the customer index)

Each txt-file is structured as follows:
- the first line contains
    - The truck dimensions (L W H) 
    - The maximum number of layers per column (S)
	- The truck capacity (K)
	- The number of items (n) 
- the remaining (n) lines all contain, for each item:
    - the item index (i)
    - the item dimensions (l_i w_i h_i) 
	- the item weight (k_i)
	- the item stackability (1 if the item is stackable, 0 if it is not)
	- the item flatness (1 if the surface of the item is flat, 0 if it is not) 
	- the customer index (c_i, only for REALISTIC_C)

Finally, the Python file ".py" can be used to visualize the optimal packing layout generated in "./PATH_AND_NAME_OUTPUT_PICTURE". One can simply run the Python code after changing the path of the picture file within the code.
