# A Simple Massive MU-MIMO Simulator for Decentralized Baseband Processing (DBP) with Feedforward Architectures
(c) 2019 Christoph Studer 
e-mail: studer@cornell.edu 

More information about our research can be found at [http://vip.ece.cornell.edu].

### Important information 

If you are using this simulator (or parts of it) for a publication, then you *must* cite our paper:

C. Jeon, K. Li, J. Cavallaro, and C. Studer, "Decentralized Equalization with Feedforward Architectures for Massive MU-MIMO. IEEE Transactions on Signal Processing, pp. 4418 -4432, Vol. 67, No. 17, July 2019

and clearly mention this in your paper.  

### How to start a simulation:

Simply run 

```sh
MIMOsim_DBP
```

which starts a simple simulation in a massive MU-MIMO system with 16 users and 256 base station antennas using 16-QAM. The simulator uses a set of predefined parameters, such as 8 antenna clusters and different algorithms (MRC, ZF, L-MMSE, feedforward partially-decentralized L-MMSE, and feedforward fully-decentralized L-MMSE). You can specify your own system, algorithm, and simulation parameters by passing your own "par" structure (see the simulator for an example). Note that we use the same parameters used to reproduce Figure 3 in the above mentioned paper; if you want to run the simulation with different parameters, then please refer to the MATLAB code for other parameter settings.
 	
We highly recommend you to execute the code step-by-step (using MATLAB's debug mode) in order to get a detailed understanding of the simulator.

### Version history
* Version 0.1 (Aug 31 2019) - studer@cornell.edu  - initial release
