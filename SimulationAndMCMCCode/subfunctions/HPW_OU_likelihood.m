function logli = HPW_OU_likelihood(Traj,parameters,z,C)
% Calculate the log-likelihood (given z, C) for a proposed hidden state move 
% i.e. log( p(Traj | parameters, z, C) )
%
% INPUT:
% Traj - the trajectory structure
% C - centre sequence
% z - hidden state sequence
% parameters - current parameter values
% prior - structure containing all the priors 
%
% OUTPUT: logli: log( p(Traj | parameters, z, C))
%
%
%
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



D=parameters(1);
D_C=parameters(2);
kappa=parameters(3);
p_esc=parameters(4);
p_trap=parameters(5);

X=Traj.Y(:,1:2);
dX=diff(X);
Dt=diff(Traj.Y(:,3));
N=length(dX);

%old looped version
% logli=zeros(N,1);
%     
% for i=1:N
%     if z(i)==0
%         logli(i)=sum(log(normpdf(dX(i,:),0,sqrt(2*D*Dt(i)))));
%     else
%         logli(i)=sum(log(normpdf(X(i+1,:),...
%             C(i,:)+(X(i,:)-C(i,:))*exp(-kappa*Dt(i)),...
%             sqrt((D/kappa)*(1-exp(-2*kappa*Dt(i)))))));
%     end  
% end


logli=zeros(1,2);
for i=1:2
    %if z=0
    Pz0=(log(normpdf(dX(:,i),0,sqrt(2*D*Dt))));
    %if z=1
    Pz1=(log(normpdf(X(2:end,i),...
        C(:,i)+(X(1:end-1,i)-C(:,i)).*exp(-kappa*Dt),...
        sqrt((D/kappa)*(1-exp(-2*kappa*Dt))))));
    
    logli(i)=sum(Pz0(z==0)) + sum(Pz1(z==1));    
end

logli=sum(logli);

end