function logPC = HPW_OU_PC(Traj,C,z,parameters,prior)
% calculate the distribution for C, given parameter and z values, 
% i.e. p(C | all parameters, z)
%
% INPUT:
% Traj - the trajectory structure
% C - centre sequence
% z - hidden state sequence
% parameters - current parameter values
% prior - structure containing all the priors 
%
% OUTPUT: logPC: log( p(C | all parameters, z))


D=parameters(1);
D_C=parameters(2);
kappa=parameters(3);
D_est=parameters(4);

Dt=diff(Traj.Y(:,3));

dC=diff(C);

N=length(C);

% logPC=zeros(N,1);
% logPC(1)=sum(log(normpdf(C(1,:),prior.mu_C,sqrt(prior.sigma_C))));
% for i=1:N-1
%     logPC(i+1)=sum(log(normpdf(dC(i,:),0,sqrt(2*Dt(i)*(D_C*z(i) + D_est*(1-z(i)))))));
% end
% logPC=sum(logPC);

logPC=zeros(1,2);
for i=1:2
    logPC(i)=sum(log(normpdf(dC(:,i),0,...
        sqrt(2*Dt(1:end-1).*(D_C*z(1:end-1) + D_est*(1-z(1:end-1)))))));
end
%add prior contribution
logPC=logPC+sum(log(normpdf(C(1,:),prior.mu_C,sqrt(prior.sigma_C))));
%sum both dimensions
logPC=sum(logPC);

end