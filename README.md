# Coupled-UZ-GRW-forspinup

<i>or, a simplified 2.5-D subsurface flow model</i>

This Matlab program combines a stand-alone 2-D phreatic aquifer model (in the following: GRW)  with a small collection of 1-D stand-alone unsaturated zone (Richards equation) models (in the following: UZ) with the aim of simulating transient groundwater tables as well as transient unsaturated pressure heads. The version presented here is setup with the aim of establishing a fast way to pre-spin-up a large scale 3-D model. 

The UZ-models ranges from the groundwater table up to the land surface. First, the UZ models are run to aquire a recharge that is then reginonlized to the full lateral model domain an used as a top boundary for the GRW-model. As a second step, the UZ model domains are resized to fit the new groundwater table, by means of changing the cell sizes of each UZ-model, while preserving the pressure profile. No iterations are  performed here, which makes the model rather inexact and approximate (and by no means mass-conservative), but very very fast. Hence, <b>USE WITH CARE! 
