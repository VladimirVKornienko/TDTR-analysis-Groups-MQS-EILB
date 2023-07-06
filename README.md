# TDTR-analysis-Groups-MQS-EILB
This repository contains the MATLAB(c) scripts for the analysis of time-domain thermoreflectance (TDTR) experimental data. The repository is created and maintained by [MQS](https://www.aalto.fi/en/department-of-electronics-and-nanoengineering/micro-and-quantum-systems) and [EILB](https://www.aalto.fi/en/department-of-electrical-engineering-and-automation/electronics-integration-and-reliability) group members of Aalto University.

The code has been originally created in the group of Prof. David G. Cahill, and is available at:\
https://cahill.matse.illinois.edu/software-and-data/

Major non-trivial changes:\
(1) In our set-up, the mechanical delay is introduced in the probe path, not in the pump path; hence, exponential factor is replaced with '1'.\
(2) We use 400 nm light for pumping and 800 nm for probing; 'AbsProfile' variable should be modified to take into account different field propagation depth for pump and probe. This needs more investigation...\
(3) Same applies to the heating model: whether it is a top surface of the transducer that has constant temperature, or whether it is a bulk of the transducer that is nbeing heated with hot electrons excited by the pump light.

Please **do not put the set-up control software here**, due to potential copyright issues with the lock-in and delay stage manufacturers.
