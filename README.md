# Coupled-UZ-GRW-forspinup

The program combines a 2-D phreatic aquifer model (in the following: GRW)  with a small collection of 1-D unsaturated zone (Richards equation) models (in the following: UZ). 

The UZ-models ranges from the groundwater table up to the land surface. First, the UZ models are run to aquire a recharge that is then reginonlized to the full lateral model domain an used as a top boundary for the GRW-model. As a second step, the UZ model domains are resized to fit the new groundwater table, by means of changing the cell sizes of each UZ-model, while preserving the pressure profile. No iterations are  performed here, which makes the model rather inexact and approximate (and by no means mass-conservative), but very very fast. Hence, <b>USE WITH CARE! 
