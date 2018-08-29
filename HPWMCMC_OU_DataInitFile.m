
%load the trajectories
load('TrajectoryData/20gold_SLBglass_0.03%GM1Processed.mat')
%load('TrajectoryData/40gold_SLBmica_0.03%GM1Processed.mat')

%number of parallel chains to run for each trajectory
n_chains=2;


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
alg_parameters.MCMC_steps=1000;
alg_parameters.burn_in=500;
%sample rate from MCMC chains
alg_parameters.thin=1;
%number of bins for parameter histograms
alg_parameters.bins=20;
%min and max block size for {z,C} MH move
alg_parameters.MinBlockSize=2;
alg_parameters.MaxBlockSize=5;
%option to initialise p_esc,p_trap from Beta(1,1) 
%(i.e. overdispered starting points for Gelman stat)
alg_parameters.overdisp=1;
alg_parameters.p_max_overdisp=0.05;
% block sizes for blocked MH move for z and C, and blocked move for C
alg_parameters.block_options.MinBlockSize=2;
alg_parameters.block_options.MaxBlockSize=1000;
alg_parameters.block_options.MultipleBlocks=1;
        
%for debugging etc. choose which MCMC moves to turn on/off
%[D, D_C, kappa, p_esc, p_trap, single z update,single C update,MH move for z and C,blocked C move]
%onchains for MCMC run used in paper: [1 1 1 1 1 0 0 1 1]
onchains=[1 1 1 1 1 0 0 1 1];

%initial values for MCMC, use [] to sample from prior
initial_values=[];




NTraj=length(Traj);
MCMCOutputSummary=cell(n_chains,1);
MCMCOutput=cell(n_chains,1);
%run HPW MCMC algorithm 


TrajIndex = 2;


   
disp(['running on ' num2str(length(TrajToRun)) ' trajectories'])
    

for j=1:n_chains
    disp(['starting chain ' num2str(j) ' on trajectory ' num2str(TrajIndex)])
    [MCMCOutput{j},MCMCOutputSummary{j}]=HPWMCMC_OU(Traj{TrajIndex},alg_parameters,prior,onchains,initial_values);
end


% plot the MCMC output
plot_HPW_OU_MCMC_output(MCMCOutput)









