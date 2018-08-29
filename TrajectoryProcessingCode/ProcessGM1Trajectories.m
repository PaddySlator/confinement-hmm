%Process raw 20nm AuNP GM1 trajectory data by removing clear artifactual
%displacements and subsampling at rate 10.
%See "additional preprocessing of trajectories" in S1 Text for full
%details.

%Paddy Slator, Warwick Systems Biology, 12/2015

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