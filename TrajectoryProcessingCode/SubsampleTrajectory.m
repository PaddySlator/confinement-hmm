function TrajSubsampled = SubsampleTrajectory(Traj,Rate)
%Function for subsampling 3-vector trajectory 
%INPUT: Traj - structure containing Y - N by 3 trajectory vector (xpos,ypos,time)
%Rate - subsampling rate
%OUTPUT: TrajSubsampled - structure containing Y, the subsampled Trajectory



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




%keep any additional fields in the input stucture Traj (e.g. hidden state,
%simulation details)
TrajSubsampled=Traj;

%if subsampling rate is 1 - no subsampling required
if Rate==1
    TrajSubsampled.Y=Traj.Y;
    return
end

%number of timepoints 
N=length(Traj.Y);

%subsample trajectory at required rate
Y=zeros(floor(N/Rate),3);
%check if hidden states are included (i.e. if Traj is a simulation)
if isfield(Traj,'z') 
   z=zeros(1,floor(N/Rate));
end
if isfield(Traj,'C')
   C=zeros(floor(N/Rate),2);
end
j=1;
for i=1:N/Rate-1
    Y(j,:)=Traj.Y(1+(i-1)*Rate,:);
    if isfield(Traj,'z')
        z(j)=Traj.z(1+(i-1)*Rate);
    end
    if isfield(Traj,'C')
       C(j,:)=Traj.C(1+(i-1)*Rate,:); 
    end
    
    j=j+1;
end
i=floor(N/Rate);
Y(j,:)=Traj.Y(1+(i-1)*Rate,:);
if isfield(Traj,'z')
    z(j)=Traj.z(1+(i-1)*Rate);
end
if isfield(Traj,'C')
    C(j,:)=Traj.C(1+(i-1)*Rate,:);
end
    
TrajSubsampled.Y=Y;
TrajSubsampled.z=z;
TrajSubsampled.C=C;

end