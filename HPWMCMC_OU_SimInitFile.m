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

% %% plot the MCMC output
% %parameter histograms and chains
% Parameters=MCMCOutputSummary{1}.Parameters;
% ParameterLabels={'D','D_C','\kappa','p_{esc}','p_{trap}'};
% for i=1:length(Parameters)
%     figure;
%     subplot(1,2,1);hold on;
%     for j=1:n_chains
%         histogram(MCMCOutputSummary{j}.ParameterChains(alg_parameters.HPW.burn_in+1:end,i))
%     end
%     plot(SimulatedTraj.parameters(i),0,'o')
%     xlabel(ParameterLabels{i})
%     ylabel('Frequency')   
%     legend_string=[];
%     for j=1:n_chains
%         legend_string{j} = ['MCMC samples: chain ' num2str(j)];
%     end
%     legend_string{end + 1} = 'simulated value';
%     %legend('MCMC samples','simulated value')
%     legend(legend_string)
%     subplot(1,2,2);hold on;
%     for j=1:n_chains
%         plot(MCMCOutputSummary{j}.ParameterChains(:,i))
%     end
%     xlabel('MCMC step')
%     ylabel(ParameterLabels{i})
% end
% 
% %hidden states
% %z
% figure; hold on;
% plot(SimulatedTraj.Y(1:end-1,3),SimulatedTraj.z,'k','LineWidth',2)
% for j=1:n_chains
%     plot(SimulatedTraj.Y(1:end-1,3),MCMCOutput{j}.z_mean,'--','LineWidth',1)
% end
% legend_string = {'simulated z'};
% for j=1:n_chains
%     legend_string{j+1} = ['mean inferred z: chain ' num2str(j)];
% end
% legend(legend_string)
% xlabel('Time (s)')
% ylabel('Confinement probability')
% 
% %C
% figure;
% subplot(1,2,1);hold on;
% plot(SimulatedTraj.Y(1:end-1,3),SimulatedTraj.C(:,1));
% for j=1:n_chains
%     plot(SimulatedTraj.Y(1:end-1,3),MCMCOutput{j}.mean_Cx);
% end
% legend_string = {'simulated C_1'};
% for j=1:n_chains
%     legend_string{j+1} = ['mean inferred C_1: chain ' num2str(j)];
% end
% legend(legend_string)
% xlabel('Time (s)')
% ylabel('C_1')
% 
% subplot(1,2,2);hold on;
% plot(SimulatedTraj.Y(1:end-1,3),SimulatedTraj.C(:,1));
% for j=1:n_chains
%     plot(SimulatedTraj.Y(1:end-1,3),MCMCOutput{j}.mean_Cx);
% end
% legend_string = {'simulated C_2'};
% for j=1:n_chains
%     legend_string{j+1} = ['mean inferred C_2: chain ' num2str(j)];
% end
% legend(legend_string)
% xlabel('Time (s)')
% ylabel('C_2')
% 
% 
% %trajectory coloured by inferred z (confinement state)
% for j=1:n_chains
%     figure;hold on;axis off;
%     X=SimulatedTraj.Y(1:end-1,1)';
%     Y=SimulatedTraj.Y(1:end-1,2)';
%     Z=zeros(size(X));
%     col=MCMCOutput{j}.z_mean;
%     caxis([0 1])
%     surface([X;X],[Y;Y],[Z;Z],[col;col],'facecol','no','edgecol','interp','linew',.5);
%     colorbar('Location','SouthOutside','Ticks',[0 1],'TickLabels',{'free','confined'})
%     title(['chain ' num2str(j)]);
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
