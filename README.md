# EKMC
ensemble kerneled matrix completion
The file ‘toy_exam_data’ contains the high-resolution data used in the toy example.
The file ‘data’ consists of synthetic data simulated by SUMO, where ‘outloop1’ and ‘outloop2’ denote the 1st and 2nd day’s data.
The m.files ‘rbf1’,’nlmc_bcd’ are basic functions shared by the toy example and the synthetic data (SUMO data) based experiments. Specifically, ‘rbf1’ is the kernel function and ‘nlmc_bcd’ is the code for block coordinate descent algorithm.
The folder ‘code-toy’ contains specific m files for the toy example, where ‘kmc’ is the code for kernelized matrix completion and ‘ekmc_num’ is the code used to the ensembled kmc algorithm for the toy example.
The folder ‘code-2’ contains files for the SUMO data based experiments. Specifically, ‘load_SUMO’ is the code used to load SUMO data and ‘ekmc_SUMO’ contains the code used to run the ensembled kmc on SUMO data. The difference between ‘kmc’ and ‘kmc1’ is they have different outputs.
The folder ‘SUMO’ contains the files for the simulation in SUMO:
Input files of dfrouter:
-dettype: defines the ids and positions of detectors;
-meaflows: defines the flows;
-network.net.xml: the network file;
Output files of dfrouter
-vehicle.xml: Saves vehicle positions as pois to FILE;
-routes.rou: Saves routes to FILE;
Input files of sumo-gui:
-vehicle.xml: - -;
-routes.rou: - -;
-network.sumo.cfg: execution file in sumo-gui;
-loopdet.xml: additional file of detector definitions;
Output files:
-outloop.xml: Saves the detector status to FILE;

Notes:
(1) All m.files should be added to the matlab path.
(2) Please refer to the https://sumo.dlr.de/docs/dfrouter.html for details on running SUMO.
