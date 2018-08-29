This repository contains data and analysis code relating to the Biophysical Journal submission: 
A Hidden Markov Model for Detecting Confinement in Single Particle Tracking Trajectories, by Slator PJ and Burroughs NJ.

Author: Paddy Slator (p.slator[at symbol]ucl.ac.uk), Centre for Medical Image Computing, University College London
[Previously Warwick Systems Biology Centre]

Date: August 2018

There are three folders: **TrajectoryData**, **SimulationAndMCMCCode** and **TrajectoryProcessingCode**;
and two example scripts: **HPWMCMCSimInitFile.m** and **HPWMCMCDataInitFile.m**.

# Example scripts

* **HPWMCMC_OU_SimInitFile.m** Simulates the harmonic potential well (HPW) model with the same parameters as the paper simulation (Figures 1-3), runs the HPW MCMC algorithm on the simulated trajectory, and plots MCMC output.

* **HPWMCMC_OU_DataInitFile.m** runs HPW MCMC algorithm on a single trajectory from the GM1 dataset, and plots MCMC output.


# TrajectoryData contains:

* The raw GM1 trajectory files:
20gold_SLBglass_0.03%GM1.mat (71 trajectories)
40gold_SLBmica_0.03%GM1.mat (18 trajectories)

* The processed GM1 trajectory files (trajectories subsampled at rate 10, with additional processing of artefacts as described in S1 Text):
20gold_SLBglass_0.03%GM1Processed.mat (71 trajectories)
40gold_SLBmica_0.03%GM1Processed.mat (18 trajectories)

For full experimental details of the trajectories see:
Spillane KM, Ortega-Arroyo J, de Wit G, Eggeling C, Ewers H, Wallace MI, Kukura P. High-Speed Single-Particle Tracking of GM1 in Model Membranes Reveals Anomalous Diffusion due to Interleaflet Coupling and Molecular Pinning. Nano Letters. 2014 Sep;14(9):5390–5397.



# SimulationAndMCMCCode contains:

* SimulateHPW.m - function that simulates the HPW model (see Results and Algorithm 2 in Supporting Information).

* HPWMCMC_OU.m - function that runs the HPW MCMC algorithm on a trajectory (see Methods and Algorithm 1 in Supporting Information).

* subfunctions - a directory containing subfunctions for the MCMC algorithm



# TrajectoryProcessingCode contains:

* ProcessGM1Trajectories.m - Function which takes raw trajectories from 20gold_SLBglass_0.03%GM1.mat file and processes as described in "additional preprocessing of trajectories" in S1 Text. This gives the trajectories included in the 20gold_SLBglass_0.03%GM1Processed.mat file as output.

* SubsampleTrajectory.m - function for subsampling trajectories.

# License and Data Ownership 

**The matlab code is released under the GNU General Public License v3.0, see LICENSE.txt.**

**The data (the contents of the TrajectoryData directory) belongs to Philipp Kukura (philipp.kukura[at symbol]chem.ox.ac.uk), and is included here with his permission.**


