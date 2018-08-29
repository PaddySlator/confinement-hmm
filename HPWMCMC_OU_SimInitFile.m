%use same parameters as in previous simulation
D=0.5;
D_C=0.01;
kappa=10000;
Dt=1/5000;
p_trap=0.002;
p_esc=0.002;
N=5000;
%parameter vector
parameters=[D D_C kappa p_esc p_trap Dt];
%D_est - diffusion coefficient of the centre when X is not confined
D_est=parameters(1);

%add noise to sim
sigma=[];

%use OU model
sim_option.OU=1;
%choose if C tracks X when z==0 (unconfined)
sim_option.CentreTracking=1;

%Simulate trajectory for HPW model
SimulatedTraj=SimulateHPW(parameters,N,D_est,sigma,sim_option);


%Prior parameters: see Methods, including "Initial values and priors" 
%parameters for gamma prior on 1/D
%a_D=0,b_D=0 for flat prior
prior.a_D=0;
prior.b_D=0;
%max/min value of D
prior.D_max=4;
prior.D_min=0;

%parameters for gamma prior on 1/D_C
%a_D_C=0,b_D_C=0 for flat prior
prior.a_D_C=0;
prior.b_D_C=0;
%ratio between D_max and D_C_max, i.e. D_C_max=D_max*D_C_ratio
prior.D_C_ratio=1/100;
%min value of D_C
prior.D_C_min=0;

%parameters for Gaussian prior on kappa, 
%use mu_kappa=0, tau_kappa=0 for flat prior
prior.mu_kappa=0;
prior.tau_kappa=0;
prior.kappa_max=20000;
prior.kappa_min=0;

%parameters for Beta priors on p_esc,p_trap 
prior.a_esc=1;
prior.b_esc=1000;
prior.a_trap=1;
prior.b_trap=1000;

%parameters for C_1
%approximate middle of the focal area
prior.mu_C=[0 0];
prior.sigma_C=1;


%alg_parameters (for MCMC)
alg_parameters.MCMC_steps=800;
alg_parameters.burn_in=400;
%sample rate from MCMC chains
alg_parameters.thin=1;
%number of bins for parameter histograms
alg_parameters.bins=20;
%min and max block size for {z,C} MH move
alg_parameters.MinBlockSize=2;
alg_parameters.MaxBlockSize=200;
%option to initialise p_esc,p_trap from Beta(1,1) 
%(i.e. overdispered starting points for Gelman stat))
alg_parameters.overdisp=0;
% block sizes for blocked MH move for z and C, and blocked move for C
alg_parameters.block_options.MinBlockSize=2;
alg_parameters.block_options.MaxBlockSize=1000;
alg_parameters.block_options.MultipleBlocks=1;

%for debugging etc. choose which MCMC moves to turn on/off
%[D, D_C, kappa, p_esc, p_trap, single z update,single C update, MH move for z and C,blocked C move]
%onchains for MCMC run used in paper: [1 1 1 1 1 0 0 1 1]
onchains=[1 1 1 1 1 0 0 1 1];

%initial values for MCMC, use [] to sample from prior
initial_values=[];

clear MCMCOutputSummary
n_chains=2;

MCMCOutputSummary=cell(n_chains,1);
MCMCOutput=cell(n_chains,1);

for i=1:n_chains
    [MCMCOutput{i},MCMCOutputSummary{i}]=HPWMCMC_OU(SimulatedTraj,alg_parameters,prior,onchains,initial_values);
end


% plot the MCMC output
plot_HPW_OU_MCMC_output(MCMCOutput)



% LICENSE
% <confinement-hmm toolbox (MCMC algorithm for detecting confinement in single particle tracking data)>
% Copyright (C) <2018>  <Paddy J. Slator, p.slator@ucl.ac.uk>
%  
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.






