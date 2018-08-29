function plot_HPW_OU_MCMC_output(MCMCOutput)
% Plot HPW model MCMC sampler output
% INPUT: MCMCOutput - structure created by HPWMCMC_OU.m


n_chains = length(MCMCOutput);
alg_parameters = MCMCOutput{1}.alg_parameters;


%parameter histograms and chains
Parameters=MCMCOutput{1}.Parameters;
ParameterLabels=MCMCOutput{1}.Parameters;


Traj = MCMCOutput{1}.Traj;
if isfield(MCMCOutput{1}.Traj,'type') && strcmp(MCMCOutput{1}.Traj.type, 'Simulation')
        simflag = 1;
else
    simflag = 0;
end


for i=1:length(Parameters)
    figure;
    subplot(1,2,1);hold on;
    for j=1:n_chains
        histogram(MCMCOutput{j}.ParameterChains(alg_parameters.burn_in+1:end,i))
    end
    if simflag
        plot(Traj.parameters(i),0,'o')
    end
    xlabel(ParameterLabels{i})
    ylabel('Frequency')
    legend_string=[];
    for j=1:n_chains
        legend_string{j} = ['MCMC samples: chain ' num2str(j)];
    end
    if simflag
        legend_string{end + 1} = 'simulated value';
    end
    legend(legend_string)
    subplot(1,2,2);hold on;
    for j=1:n_chains
        plot(MCMCOutput{j}.ParameterChains(:,i))
    end
    xlabel('MCMC step')
    ylabel(ParameterLabels{i})
end

%hidden states
%z
figure; hold on;
if simflag
plot(Traj.Y(1:end-1,3),Traj.z,'k','LineWidth',2)
end
for j=1:n_chains
    plot(Traj.Y(1:end-1,3),MCMCOutput{j}.z_mean,'--','LineWidth',1)
end
legend_string=[];
for j=1:n_chains
    legend_string{j} = ['mean inferred z: chain ' num2str(j)];
end
if simflag
    legend_string(2:end+1) = legend_string;
    legend_string{1} = 'simulated z';
end
legend(legend_string)
xlabel('Time (s)')
ylabel('Confinement probability')

%C
figure;
subplot(1,2,1);hold on;
if simflag
plot(Traj.Y(1:end-1,3),Traj.C(:,1));
end
for j=1:n_chains
    plot(Traj.Y(1:end-1,3),MCMCOutput{j}.mean_Cx);
end
legend_string=[];
for j=1:n_chains
    legend_string{j} = ['mean inferred C_1: chain ' num2str(j)];
end
if simflag
    legend_string(2:end+1) = legend_string;
    legend_string{1} = 'simulated C_1';
end
legend(legend_string)
xlabel('Time (s)')
ylabel('C_1')

subplot(1,2,2);hold on;
if simflag
    plot(Traj.Y(1:end-1,3),Traj.C(:,1));
end
for j=1:n_chains
    plot(Traj.Y(1:end-1,3),MCMCOutput{j}.mean_Cx);
end
legend_string=[];
for j=1:n_chains
    legend_string{j} = ['mean inferred C_2: chain ' num2str(j)];
end
if simflag
    legend_string(2:end+1) = legend_string;
    legend_string{1} = 'simulated C_2';
end

legend(legend_string)
xlabel('Time (s)')
ylabel('C_2')


%trajectory coloured by inferred z (confinement state)
for j=1:n_chains
    figure;hold on;axis off;
    X=Traj.Y(1:end-1,1)';
    Y=Traj.Y(1:end-1,2)';
    Z=zeros(size(X));
    col=MCMCOutput{j}.z_mean;
    caxis([0 1])
    surface([X;X],[Y;Y],[Z;Z],[col;col],'facecol','no','edgecol','interp','linew',.5);
    colorbar('Location','SouthOutside','Ticks',[0 1],'TickLabels',{'free','confined'})
    title(['chain ' num2str(j)]);
end












