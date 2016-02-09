*SNS IKL*
Version 0.1.0 beta

The Inverse Kinematic Library (*SNS IKL*) contains tools to invert the differential kinematic of a robot (this version is tested only on fixed base manipulator). 
The methods implemented on the IKL are:
- STD 		the standard pseudoinversion of the jacobian
- SCALE 	pseudoinversion with task scaling to satisfy joint constraints 
- SNS		the SNS algorithm which permits to consider hard constraints at joint level

All this methods permits to solve the inversion of the kinematic of a prioritized stack of task. (note that only the SNS algorithm guarantees to satisfy all joint constraints).   


Author: Fabrizio Flacco
Dipartimento di Ingegneria Informatica, Automatica e Gestionale (DIAG)
Università di Roma "La Sapienza"
Rome, Italy
