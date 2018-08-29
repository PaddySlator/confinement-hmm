function logPz = HPW_OU_Pz(z,parameters)
% calculate the distribution for z, given parameter values, (doesn't depend on C) 
% i.e. p(z | all parameters)
%
% INPUT:
% z - hidden state sequence
% parameters - current parameter values
%
% OUTPUT: logPz: log( p(Z | all parameters))
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

p_esc=parameters(1);
p_trap=parameters(2);

N=length(z);


logPz=sum(log(...
    z(2:end).*z(1:end-1).*(1-p_esc)...
    +z(2:end).*(1-z(1:end-1)).*p_trap...
    +(1-z(2:end)).*z(1:end-1).*p_esc...
    +(1-z(2:end)).*(1-z(1:end-1)).*(1-p_trap)...
    ));
            
%add contribution from prior (based on stationary probabilities)
logPz=logPz+log(z(1)*p_trap/(p_trap+p_esc) + (1-z(1))*p_esc/(p_trap+p_esc));


end