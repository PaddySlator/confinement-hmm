function [MCMCOutput,MCMCOutputSummary] = HPWMCMC_OU(Traj,alg_parameters,prior,onchains,initial_values)
% MCMC sampler to infer model parameters and hidden states for 2D Brownian motion with a switching harmonic potential well
% using Ornstein-Uhlenbeck solution to the SDE
%
% INPUTS:
% Traj - the trajectory structure to fit to, must contain Traj.Y, the 3-vector of the trajectory (Xposition,Yposition,time) 
%
% alg_parameters structure with fields choosing MCMC algorithm options
% alg_parameters.MCMC_steps - total number of MCMC samples
% alg_parameters.burn_in - number of samples to discard
% alg_parameters.thin - sample rate for MCMC samples
% alg_parameters.bins - number of bins for parameter histograms
% alg_parameters.MinBlockSize, alg_parameters.MaxBlockSize - min and max block size for {z,C} MH move
% alg_parameters.overdisp - option to initialise p_esc,p_trap from Beta(1,1) - to give
% overdispersed starting points for Gelman stat
% alg_parameters.p_max_overdisp - option to scale the overdispered Beta(1,1) samples to a
% maximum value
% alg_parameters.block_options.MinBlockSize - minimum block size for
% blocked MH move for z and C, and blocked move for C
% alg_parameters.block_options.MaxBlockSize - maximum block size for
% blocked moves
% alg_parameters.block_options.MultipleBlocks - number of blocks to sample at each MCMC step

% prior - structure defining the priors to use 
%
% prior.a_D=0, prior.b_D=0 - parameters for gamma prior on 1/D
% (use a_D=0,b_D=0 for flat prior)
%
% prior.D_max, prior.D_min - max/min value of D
%
% prior.a_D_C, prior.b_D_C - parameters for gamma prior on 1/D_C
% (use a_D_C=0,b_D_C=0 for flat prior)
%
% prior.D_C_ratio - ratio between D_max and D_C_max, i.e. D_C_max=D_max*D_C_ratio
% prior.D_C_min=0 - min value of D_C
%
% prior.mu_kappa,prior.tau_kappa -  Gaussian prior on kappa, 
% (use mu_kappa=0, tau_kappa=0 for flat prior)
% prior.kappa_max, prior.kappa_min - truncate the prior
%
% prior.a_esc, prior.b_esc - Beta prior on p_esc
% prior.a_trap, prior.b_trap - Beta prior on p_trap
%
% prior.mu_C, prior.sigma_C - normal prior on C_1
%
% onchains - vector of flags choosing which MCMC moves to turn on/off - [onD onD_C onkappa onp_esc onp_trap onz onC]
%
%
% initial_values - optional - define the initial value of any parameter or hidden state 
% if blank initial values are sampled from the priors
%
%
% OUTPUT: MCMCOutput and MCMCOutputSummary - structures containg the detailed MCMC output


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






%output the runtime
tic;

%trajectory length
N=length(Traj.Y)-1;
%2-vector of positions
X=[Traj.Y(1:end,1), Traj.Y(1:end,2)];
X=[Traj.Y(:,1) Traj.Y(:,2)];
%2-vector of displacements
dX=diff([Traj.Y(:,1),Traj.Y(:,2)]);
%vector of timesteps
Dt=diff(Traj.Y(:,3));

MCMC_steps=alg_parameters.MCMC_steps;
burn_in=alg_parameters.burn_in;
thin=alg_parameters.thin;

%debugging
%parameters
onD=onchains(1);
onD_C=onchains(2);

onkappa=onchains(3);
onp_esc=onchains(4);
onp_trap=onchains(5);
%single updates for hidden states
onz=onchains(6);
onC=onchains(7);
onzC_MH=onchains(8);
onBlockedC=onchains(9);

if onzC_MH
    Moves=zeros(MCMC_steps,1);
    AttemptedMoves=zeros(MCMC_steps,1);
    BlockSizeStep=zeros(MCMC_steps,1);
end

%set D_est equal to maximum likelihood estimate for D
D_est=0.25*mean(sum(dX.^2,2)./Dt);
MCMCOutput.D_est=D_est;
MCMCOutputSummary.D_est=D_est;


%Output which MCMC moves are being use, and which priors are being used
NMCMCMoves=7;

chains=cell(NMCMCMoves,2);
chains{1,1}='D';
chains{2,1}='D_C';
chains{3,1}='kappa';
chains{4,1}='p_esc';
chains{5,1}='p_trap';
chains{6,1}='z';
chains{7,1}='C';

for i=1:NMCMCMoves
    chains{i,2}=onchains(i);
end
MCMCOutput.chains=chains;
MCMCOutputSummary.chains=chains;


%Initialise chains
%parameters
D_chain=zeros(MCMC_steps/thin,1);
D_C_chain=zeros(MCMC_steps/thin,1);
kappa_chain=zeros(MCMC_steps/thin,1);
p_trap_chain=zeros(MCMC_steps/thin,1);
p_esc_chain=zeros(MCMC_steps/thin,1);

Parameters={'D','D_C','kappa','p_esc','p_trap'};

%all parameters together
ParameterChains=zeros(MCMC_steps/thin,length(Parameters));

%extract prior parameters

D_min=prior.D_min;
D_max=prior.D_max;

D_C_min=prior.D_C_min;
D_C_max=prior.D_C_ratio*D_max;

MCMCOutput.D_max=D_max;

kappa_max=prior.kappa_max;
kappa_min=prior.kappa_min;


if onD
    if isfield(initial_values,'D')
        D = initial_values.D;
    else
        %sample from prior
        if prior.a_D==0 && prior.b_D==0
            %uniform
            D=D_min+rand*(D_max-D_min);
        else
            %inverse gamma
            D=1/gamrnd(prior.a_D,1/prior.b_D);
        end      
    end
else
    %for debugging
    D=Traj.parameters(1);
end

if onD_C
    if isfield(initial_values,'D_C')
        D_C = initial_values.D_C;
    else
        if prior.a_D_C == 0 || prior.b_D_C == 0
            %uniform
            D_C=D_C_min+rand*(D_C_max - D_C_min);
        else
            %inverse gamma
            D_C=1/gamrnd(prior.a_D_C,1/prior.b_D_C);
        end               
    end
else
    %for debugging
    D_C=Traj.parameters(2);
end


if onkappa
    if isfield(initial_values,'kappa')
        kappa=initial_values.kappa;
    else
        if prior.mu_kappa == 0 || prior.tau_kappa == 0
            %uniform
            kappa=kappa_min+rand*(kappa_max-kappa_min);
        else
            %normal
            kappa=prior.mu_kappa + sqrt(1/(2*prior.tau_kappa))*randn;
            while kappa < kappa_min || kappa > kappa_max
                kappa=prior.mu_kappa + sqrt(1/(2*prior.tau_kappa))*randn;
            end
        end        
    end
    %count the number of MH moves
    kappa_moves=0;
else
    %for debugging
    kappa=Traj.parameters(3);
end

%p_esc
if onp_esc
    if isfield(initial_values,'p_esc')
        p_esc=initial_values.p_esc;
    elseif (prior.a_esc || prior.b_esc) == 0
        p_esc=rand;        
    else
        p_esc = betarnd(prior.a_esc,prior.b_esc);      
        if alg_parameters.overdisp
            p_esc=rand*alg_parameters.p_max_overdisp;
        end
    end
else
    p_esc=Traj.parameters(4);
end
%p_trap
if onp_trap
    if isfield(initial_values,'p_trap')
        p_trap=initial_values.p_trap;
    elseif (prior.a_trap || prior.b_trap) == 0
        p_trap=rand;       
    else
        p_trap = betarnd(prior.a_trap,prior.b_trap);
        if alg_parameters.overdisp
            p_trap=rand*alg_parameters.p_max_overdisp;
        end
    end
else
    p_trap=Traj.parameters(5);
end

D_chain(1)=D;
D_C_chain(1)=D_C;
kappa_chain(1)=kappa;
p_trap_chain(1)=p_trap;
p_esc_chain(1)=p_esc;

MCMCOutput.initial_values.D=D;
MCMCOutput.initial_values.D_C=D_C;
MCMCOutput.initial_values.kappa=kappa;
MCMCOutput.initial_values.p_esc=p_esc;
MCMCOutput.initial_values.p_trap=p_trap;

%hidden states
z_chain=zeros(MCMC_steps/thin,N);
Cx_chain=zeros(MCMC_steps/thin,N);
Cy_chain=zeros(MCMC_steps/thin,N);

% confinement state
if onz || onzC_MH
    if isfield(initial_values,'z')
        z=initial_values.z;
    else
        z=zeros(N,1);
        %flat prior on confinement of first point
        z(1)=randsample([0 1],1);
        for i=2:N
            if z(i-1)==1 && rand > p_esc
                z(i)=1;
            end
            if z(i-1)==0 && rand < p_trap
                z(i)=1;
            end
        end
    end
else
    z=Traj.z';
end
z_chain(1,:)=z;



%C
if onC || onzC_MH || onBlockedC
    if isfield(initial_values,'C')
        C=initial_values.C;
        if ~(onU || onU_block)
            X=X;
        end
    else %initialise C by smoothing X
        C(:,1)=smooth(X(1:end-1,1))';
        C(:,2)=smooth(X(1:end-1,2))';
        disp('Initialised C by smoothing X')
    end   
else
   C=Traj.C; 
end
Cx_chain(1,:)=C(:,1);
Cy_chain(1,:)=C(:,2);

MCMCOutput.initial_values.z=z;
MCMCOutput.initial_values.C=C;
    


for step=2:MCMC_steps
    %Update parameters
    
    %Gibbs update for D
    if onD
        D = 1/(gamrnd(prior.a_D + N - 1,...
            1/(prior.b_D + ...
            0.5*sum(sum(( dX - [z z].*(C-X(1:end-1,:)+(X(1:end-1,:)-C).*exp([-kappa*Dt -kappa*Dt]))).^2,2)...
            ./(2*(1-z).*Dt + (z./kappa).*(1-exp(-2*kappa*Dt)))))));
                        
        
        %enforce max and min D condition
        if D > D_max || D < D_min
            D=D_chain(step-1);
        end
    end
    
    %Gibbs update for D_C
    if onD_C
        %2-vector of centre displacements
        dC=diff(C);
        if sum(z)==0
            %if there are no confined steps (sum(z)==0) then sample from the prior
            D_C=D_C_min+rand*(D_C_max-D_C_min);
        else
            D_C = 1/(gamrnd(prior.a_D_C + sum(z) + 1 , ...
                1/(prior.b_D_C + 0.25*sum(z(1:end-1).*sum(dC.^2,2)./Dt(1:end-1)))));
        end
        %enforce max and min D_C condition
        if D_C > D_C_max || D_C < D_C_min
            D_C=D_C_chain(step-1);
        end
    end
        
    %kappa update
    if onkappa
        S_kappa=400;        
        kappa_prop=kappa+S_kappa*randn;
                           
       % calculate density proportional to  \pi(kappa|...)
        p_kappa=sum(log(normpdf(dX,...
            C-X(1:end-1,:)+(X(1:end-1,:)-C).*exp([-kappa*Dt -kappa*Dt]),...
            repmat(sqrt(D*(1-exp(-2*kappa*Dt(i)))/kappa),N,2))),2);
                        
        p_kappa_prop=sum(log(normpdf(dX,...
            C-X(1:end-1,:)+(X(1:end-1,:)-C).*exp([-kappa_prop*Dt -kappa_prop*Dt]),...
            repmat(sqrt(D*(1-exp(-2*kappa_prop*Dt(i)))/kappa_prop),N,2))),2);
                                                
        %set p_kappa and p_kappa_prop to 0 where z=0 
        %(can't think of an easy way of vectorising without writing out the whole normal pdf)        
        p_kappa(z==0)=0;
        p_kappa_prop(z==0)=0;    
                      
        if (kappa_prop > kappa_min && kappa_prop < kappa_max) %enforce max and min kappa condition
            if sum(p_kappa_prop)>sum(p_kappa) || log(rand)<(sum(p_kappa_prop)-sum(p_kappa)) %MH acceptance rate
                kappa=kappa_prop;  
                kappa_moves=kappa_moves+1;
            end
        end
    end
    
    if onp_esc || onp_trap
        %number of trapping and escape events
        n_11=sum(z(1:end-1)==1 & z(2:end)==1);
        n_00=sum(z(1:end-1)==0 & z(2:end)==0);
        n_01=sum(z(1:end-1)==0 & z(2:end)==1);
        n_10=sum(z(1:end-1)==1 & z(2:end)==0);
    end
    %p_esc
    if onp_esc
        p_esc=betarnd(prior.a_esc + n_10  , prior.b_esc + n_11);
    end
    %p_trap
    if onp_trap
        p_trap=betarnd(prior.a_trap + n_01, prior.b_trap + n_00);
    end
        
    %z update (confinement state)
    if onz
        %dX=diff(X);
        dC=diff(C);
        pi_trap=p_trap/(p_trap+p_esc);
        pi_esc=1-pi_trap;
        i=1;
        
        p_z0 = pi_esc*...
            prod(normpdf(dX(i,:),0,sqrt(2*D*Dt(i))))...
            *prod(normpdf(dC(i,:),0,sqrt(2*D_est*Dt(i))))... 
            *(z(i+1)*p_trap + (1-z(i+1))*(1-p_trap));
        
        p_z1 = pi_trap*...
            prod(normpdf(dX(i,:),...
                C(i,:)-X(i,:)+(X(i,:)-C(i,:))*exp(-kappa*Dt(i)),...
                sqrt(D*(...
                    (1-z(i))*2*Dt(i)...
                    +(z(i)/kappa)*(1-exp(-2*kappa*Dt(i)))...
                ))))...
            *prod(normpdf(dC(i,:),0,sqrt(2*D_C*Dt(i))))... 
            *(z(i+1)*(1-p_esc) + (1-z(i+1))*p_esc);
        
        if rand < p_z1/(p_z1+p_z0)
            z(i)=1;
        else
            z(i)=0;
            if sum(z) <= 2
                z(i)=1;
            end
        end
        
        PointsToUpdate=2:N-1;
        
        for i=PointsToUpdate
            
            p_z1=(z(i-1)*(1-p_esc)+ (1-z(i-1))*p_trap)...
                *prod(normpdf(dX(i,:),...
                C(i,:)-X(i,:)+(X(i,:)-C(i,:))*exp(-kappa*Dt(i)),...
                sqrt(D*(...
                (1-z(i))*2*Dt(i)...
                +(z(i)/kappa)*(1-exp(-2*kappa*Dt(i)))...
                ))))...
                *prod(normpdf(dC(i,:),0,sqrt(2*D_C*Dt(i))))...
                *(z(i+1)*(1-p_esc)+ (1-z(i+1))*p_esc);
                       
%             mu_z1=C(i,:)-X(i,:)+(X(i,:)-C(i,:))*exp(-kappa*Dt(i));
%             sigma_z1=sqrt(D*(...
%                 (1-z(i))*2*Dt(i)...
%                 +(z(i)/kappa)*(1-exp(-2*kappa*Dt(i)))));
%             
%             p_z1=(z(i-1)*(1-p_esc)+ (1-z(i-1))*p_trap)...
%                 *prod(exp(-0.5 * ((dX(i,:) - mu_z1)./sigma_z1).^2) ./ (sqrt(2*pi) .* sigma_z1))...
%                 *prod(exp(-0.5 * (dC(i,:)./sqrt(2*D_C*Dt(i))).^2) ./ (sqrt(2*pi) .* sqrt(2*D_C*Dt(i))))...              
%                 *(z(i+1)*(1-p_esc)+ (1-z(i+1))*p_esc);
                         
            p_z0=(z(i-1)*p_esc + (1-z(i-1))*(1-p_trap))...
                *prod(normpdf(dX(i,:),0,sqrt(2*D*Dt(i))))...
                *prod(normpdf(dC(i,:),0,sqrt(2*D_est*Dt(i))))...  
                *(z(i+1)*p_trap + (1-z(i+1))*(1-p_trap));
                                    
%             p_z0=(z(i-1)*p_esc + (1-z(i-1))*(1-p_trap))...
%                 *prod(exp(-0.5 * (dX(i,:)./sqrt(2*D*Dt(i))).^2) ./ (sqrt(2*pi) .* sqrt(2*D*Dt(i))))...                
%                 *prod(exp(-0.5 * (dC(i,:)./sqrt(2*D_est*Dt(i))).^2) ./ (sqrt(2*pi) .* sqrt(2*D_est*Dt(i))))...                 
%                 *(z(i+1)*p_trap + (1-z(i+1))*(1-p_trap));
                                    
            if rand < p_z1/(p_z1+p_z0)
                z(i)=1;
            else
                z(i)=0;
                if sum(z) <= 2
                    z(i)=1;
                end
            end
        end
        i=N;
        
        p_z0 = (z(i-1)*p_esc + (1-z(i-1))*(1-p_trap))...
            *prod(normpdf(dX(i,:),0,sqrt(2*D*Dt(i))));
        
        p_z1=(z(i-1)*(1-p_esc)+ (1-z(i-1))*p_trap)...
            *prod(normpdf(dX(i,:),...
                    C(i,:)-X(i,:)+(X(i,:)-C(i,:))*exp(-kappa*Dt(i)),...
                    sqrt(D*(...
                        (1-z(i))*2*Dt(i)...
                        +(z(i)/kappa)*(1-exp(-2*kappa*Dt(i)))...
                    ))));
        
        
        if rand < p_z1/(p_z1+p_z0)
            z(i)=1;
        else
            z(i)=0;
            if sum(z) <= 2
                z(i)=1;
            end
        end
    end
    
    %C update
    if onC
        %dX=diff(X);
        tau_C=zeros(N,1);
        mu_C=zeros(N,2);
        
        i=1;
        tau1=(z(i)*(1-exp(-kappa*Dt(i)))^2)/...
            (D*(1-exp(-2*kappa*Dt(i)))/kappa);
        tau3=1/(2*Dt(i)*(D_C*z(i)+D_est*(1-z(i))));
        
        mu1=z(i)*(dX(i,:)+X(i,:)*(1-exp(-kappa*Dt(i)))/...
            (1-exp(-kappa*Dt(i))));
        mu3=C(i+1,:);
        
        tau_C(i)=tau1+tau3;
        mu_C(i,:)=(mu1*tau1 + mu3*tau3)*(1/tau_C(i));
        
        C(i,:)= sqrt(1/tau_C(i))*randn(1,2)...
            +mu_C(i,:);
        
        for i=2:N-1
            tau1=(z(i)*(1-exp(-kappa*Dt(i)))^2)/...
                (D*(1-exp(-2*kappa*Dt(i)))/kappa);
            tau2=1/(2*Dt(i-1)*(D_C*z(i-1)+D_est*(1-z(i-1))));
            tau3=1/(2*Dt(i)*(D_C*z(i)+D_est*(1-z(i))));
            
            mu1=z(i)*(dX(i,:)+X(i,:)*(1-exp(-kappa*Dt(i)))/...
                (1-exp(-kappa*Dt(i))));
            mu2=C(i-1,:);
            mu3=C(i+1,:);
            
            tau_C(i)=tau1+tau2+tau3;
            
            mu_C(i,:)=(mu1*tau1 + mu2*tau2 + mu3*tau3)*(1/tau_C(i));
            
            C(i,:)= sqrt(1/tau_C(i))*randn(1,2)...
                +mu_C(i,:);
        end
        i=N;
        
        tau1=(z(i)*(1-exp(-kappa*Dt(i)))^2)/...
            (D*(1-exp(-2*kappa*Dt(i)))/kappa);
        tau2=1/(2*Dt(i-1)*(D_C*z(i-1)+D_est*(1-z(i-1))));
        
        mu1=z(i)*(dX(i,:)+X(i,:)*(1-exp(-kappa*Dt(i)))/...
            (1-exp(-kappa*Dt(i))));
        mu2=C(i-1,:);
        
        tau_C(i)=tau1+tau2;
        mu_C(i,:)=(mu1*tau1 + mu2*tau2)*(1/tau_C(i));
        C(i,:)= sqrt(1/tau_C(i))*randn(1,2)...
            +mu_C(i,:);                
    end
           
    %MH move for z and C
    if onzC_MH        
        %sample the block size and blocks to move        
        Blocks=SampleBlockToMove(N,alg_parameters.block_options); 
        
        for k=1:size(Blocks,1)
            %occasionally get a CPrecisionMatrix which is close to singular -
            %so either get a warning or sometimes an error
            %only seems to happen if there is a single z==1 (i.e. one timepoint is confined)
            %if this happens just skip the z,C update (doesn't affect posterior
            %as we are just rearraging the update order)
            try
                options.BlockToMove=Blocks(k,1):Blocks(k,2);
                options.ProposalOptions='Blocked';
                %propose MH move (on this block) for z, and calculate p(z|z') and p(z'|z)
                [zProp,LogzDensity,LogzPropDensity]=...
                    MetropolisHastingsForHiddenState(z,options);
                
                %propose value of C (changing C only on this block) using blocked Gibbs move
                options.sample=1;
                [CProp,LogCPropDensity]=HPWMCMC_OU_BlockedCMove(Traj,C,zProp,[D D_C kappa D_est],prior,options);
                %calculate the Gibbs density of the reverse move (CProp | ZProp -> C)
                options.sample=0;
                options.CForDensity=C;
                %(doesn't actually make any different if you use C or CProp in this function
                %- only need the two timepoints bordering the block, which are the same in both.)
                [~,LogCDensity]=HPWMCMC_OU_BlockedCMove(Traj,CProp,z,[D D_C kappa D_est],prior,options);
                
                %calculate p(z,C|X,theta)
                %\pi(C|z,theta,X)
                PC=HPW_OU_PC(Traj,C,z,[D D_C kappa D_est],prior);
                PropPC=HPW_OU_PC(Traj,CProp,zProp,[D D_C kappa D_est],prior);
                %\pi(z|theta)
                Pz=HPW_OU_Pz(z,[p_esc p_trap]);
                PropPz=HPW_OU_Pz(zProp,[p_esc p_trap]);
                %\pi(X|z,C,theta)
                li=HPW_OU_likelihood(Traj,[D D_C kappa p_esc p_trap],z,C);
                Propli=HPW_OU_likelihood(Traj,[D D_C kappa p_esc p_trap],zProp,CProp);
                
                Paccept=min(log(1),...
                    (Propli+PropPz+PropPC+LogzDensity+LogCDensity)...
                    -(li+Pz+PC+LogzPropDensity+LogCPropDensity));
                
                if log(rand) < Paccept
                    Moves(step)=Moves(step)+1;
                    
                    z=zProp;
                    C=CProp;
                end
            catch
                
            end
        end
    end
    
    if onBlockedC      
        %sample the block size and blocks to move        
        Blocks=SampleBlockToMove(N,alg_parameters.block_options); 
        
        for k=1:size(Blocks,1)
            options.sample=1;
            options.BlockToMove=Blocks(k,1):Blocks(k,2);
            %occasionally get a CPrecisionMatrix which is close to singular -
            %so either get a warning or sometimes an error
            %only seems to happen if there is a single z==1 (i.e. one timepoint is confined)
            %if this happens just skip the C update (doesn't affect posterior
            %as we are just rearraging the update order)
            try
                C=HPWMCMC_OU_BlockedCMove(Traj,C,z,[D D_C kappa D_est],prior,options);
            catch
                
            end
        end
    end
    
    
    %Update chains, sample at the sampling rate
    if mod(step,thin) == 0
        %parameters
        D_chain(step/thin)=D;
        D_C_chain(step/thin) = D_C;
        kappa_chain(step/thin) = kappa;
        p_trap_chain(step/thin)=p_trap;
        p_esc_chain(step/thin)=p_esc;
        
        ParameterChains(step/thin,:)=[D D_C kappa p_trap p_esc];
        
        %hidden states
        z_chain(step/thin,:)=z;
        Cx_chain(step/thin,:) = C(:,1);
        Cy_chain(step/thin,:) = C(:,2);                        
    end                 
end


%Save the final state (can run another chain using these as initial values)
%parameters
final_state.D=D;
final_state.D_C=D_C;
final_state.kappa=kappa;
final_state.p_esc=p_esc;
final_state.p_trap=p_trap;
%hidden states
final_state.Z=X;
final_state.C=C;
final_state.z=z;


%Two output structures - MCMCOutput and MCMCOutputSummary

%MCMCOutputSummary (smaller structure)

%parameter chains
MCMCOutputSummary.D_chain=D_chain;
MCMCOutputSummary.D_C_chain=D_C_chain;
MCMCOutputSummary.kappa_chain=kappa_chain;
MCMCOutputSummary.p_trap_chain=p_trap_chain;
MCMCOutputSummary.p_esc_chain=p_esc_chain;
MCMCOutputSummary.ParameterChains=ParameterChains;
    
%samples from parameter posteriors
MCMCOutputSummary.D_PosteriorSamples=D_chain((burn_in/thin):end);
MCMCOutputSummary.D_C_PosteriorSamples=D_C_chain((burn_in/thin):end);
MCMCOutputSummary.kappa_PosteriorSamples=kappa_chain((burn_in/thin):end);
MCMCOutputSummary.p_trap_PosteriorSamples=p_trap_chain((burn_in/thin):end);
MCMCOutputSummary.p_esc_PosteriorSamples=p_esc_chain((burn_in/thin):end);

MCMCOutputSummary.ParameterPosteriorSamples=[MCMCOutputSummary.D_PosteriorSamples...
    MCMCOutputSummary.D_C_PosteriorSamples...
    MCMCOutputSummary.kappa_PosteriorSamples...
    MCMCOutputSummary.p_esc_PosteriorSamples...
    MCMCOutputSummary.p_trap_PosteriorSamples];


%mean hidden states
%C
MCMCOutputSummary.mean_Cx=mean(Cx_chain(burn_in/thin:end,:));
MCMCOutputSummary.mean_Cy=mean(Cy_chain(burn_in/thin:end,:));
MCMCOutputSummary.sd_Cx=std(Cx_chain(burn_in/thin:end,:));
MCMCOutputSummary.sd_Cy=std(Cy_chain(burn_in/thin:end,:));
%z
MCMCOutputSummary.z_mean=mean(z_chain(burn_in/thin:end,:));
MCMCOutputSummary.z_sd=std(z_chain(burn_in/thin:end,:));



%Simulation
MCMCOutputSummary.Traj=Traj;
%Algorithm parameters
MCMCOutputSummary.alg_parameters=alg_parameters;
%prior
MCMCOutputSummary.prior=prior;
%Final parameter and hidden state values (useful for carrying on chains)
MCMCOutputSummary.final_state=final_state;

MCMCOutputSummary.Model='MCMC run for HPW model';
MCMCOutputSummary.Parameters={'D','D_C','kappa','p_esc','p_trap'};

%Metropolis-Hastings stats
MCMCOutputSummary.Moves=Moves;
MCMCOutputSummary.kappa_moves=kappa_moves;

%MCMCOutput (much bigger structure)
MCMCOutput=MCMCOutputSummary;

%add hidden state chains (which take up a lot of memory!)
MCMCOutput.Cx_chain=Cx_chain;
MCMCOutput.Cy_chain=Cy_chain;
MCMCOutput.z_chain=z_chain;



toc;
end