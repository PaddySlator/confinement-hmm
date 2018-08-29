%Process raw 20nm AuNP GM1 trajectory data by removing clear artifactual
%displacements and subsampling at rate 10.
%See "additional preprocessing of trajectories" in S1 Text for full
%details.

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


%load Traj - cell containing all trajectories
load('TrajectoryData/20gold_SLBglass_0.03%GM1.mat')

NTraj=length(Traj);

%truncate trajectories with clear artifactual displacements 
TrajToTruncate=[9 13 29 33 35 51 52];

%trajectories with artifacts towards beginning of trajectory
KeepEnd=9;
%trajectories with artifacts towards end of trajectory
KeepStart=[13 29 33 35 51 52];

%point where artifacts either begin or end
TruncationPoint=[1151 16570 40369 46130 47249 10155 39795];

TrajTruncated=Traj;
%Truncate trajectories 
for i=1:length(TrajToTruncate)
    if any(TrajToTruncate(i)==KeepEnd)
        %keep data after truncation point
        TrajTruncated{TrajToTruncate(i)}.Y=...
            Traj{TrajToTruncate(i)}.Y(TruncationPoint(i):end,:);
    end
    if any(TrajToTruncate(i)==KeepStart)
        %keep data before truncation point
        TrajTruncated{TrajToTruncate(i)}.Y=...
            Traj{TrajToTruncate(i)}.Y(1:TruncationPoint(i),:);
    end
end

%Subsample trajectories at rate 10 (S1 Text for explanation)
SubsamplingRate=10;
TrajSubsampled=cell(1,NTraj);
for i=1:NTraj
    TrajSubsampled{i}=SubsampleTrajectory(TrajTruncated{i},SubsamplingRate);
end