function Traj = SimulateHPW(parameters,N,D_est,sigma,option)
%Simulate 2D Brownian motion with switching harmonic well (HPW model)
%INPUT: parameters - [D D_C kappa p_esc p_trap Dt],
% N - number of timesteps,
% D_est - diffusion rate for centre when X unconfined 
% sigma - standard deviation of static measurement error, use [] for
% noiseless simulation

% option - various simulation options  
% option.OU - SDE approximation flag: 1 = Ornstein-Uhlenbeck, 0 = Euler
% option.fixed_z - optional user defined confinement state sequence
% option.CentreTracking - flag, if =1 then C tracks X when z==0
%

%OUTPUT: Traj - Structure containing:
% Y - 3-vector of trajectory (Xposition,Yposition,time)
% z - confinement state sequence 
% C - harmonic potential well centre sequence

%Paddy Slator, Warwick Systems Biology Centre, 11/2015

%default to Ornstein-Uhlenbeck solution of SDEs
if ~isfield(option,'OU')
    option.OU = 1;
end


D=parameters(1);
D_C=parameters(2);
kappa=parameters(3);
p_esc=parameters(4);
p_trap=parameters(5);

%if no D_est specified, set to D.
if isempty(D_est)
    D_est=D;
end 
Traj.D_est=D_est;

% Equal time steps
Dt=parameters(6)*ones(N,1);


%simulate confinement state z
if isfield(option,'fixed_z') && ~isempty(option.fixed_z)
    z=option.fixed_z;
else
    %(Markov chain with transition probabilities p_trap,p_esc)
    z=zeros(N,1);
    %simulate z(1) using the stationary probability
    if rand < p_trap/(p_trap+p_esc)
        z(1)=1;
    end
    for i=2:N
        if z(i-1)==1 && rand > p_esc
            z(i)=1;
        end
        if z(i-1)==0 && rand < p_trap
            z(i)=1;
        end
    end
end

%centre vector
C=zeros(N,2);
%particle position vector
X=zeros(N+1,2);
%X,C start at (0,0)
i=1;
C(i,:)=[0 0];
X(i,:)=[0 0];

for i=1:N-1
    if option.CentreTracking        
        %C tracks X when z==0
        %C(i+1,:)= C(i,:)...
        %    + (1-z(i))*kappa*Dt(i)*(sum(X)-C(i,:))...
        %    + sqrt(2*Dt(i)*(z(i)*D_C...
        %    + (1-z(i))*D_est))*randn(1,2);
        if z(i)==1
            C(i+1,:)=C(i,:)+sqrt(2*Dt(i)*D_C)*randn(1,2);
        else
            %just keep C close to X
            C(i+1,:)=X(i,:)+sqrt(2*Dt(i)*D_C)*randn(1,2);            
        end
    else              
        %C diffuses when z==0
        %Euler approximation
        C(i+1,:)= C(i,:)...
            + sqrt(2*Dt(i)*(z(i)*D_C...
            + (1-z(i))*D_est))*randn(1,2);
    end
    
    if option.OU
        %use Ornstein-Uhlenbeck solution
        if z(i)==1
            X(i+1,:)= C(i,:)...
                + (X(i,:)-C(i,:))*exp(-kappa*Dt(i))...
                + sqrt((D/kappa)*(1-exp(-2*kappa*Dt(i))))*randn(1,2);
        else
            X(i+1,:)=X(i,:)+sqrt(2*D*Dt(i))*randn(1,2);
        end
    else
        %Euler approximation to HPW model SDEs
        X(i+1,:)= X(i,:)...
            -(z(i))*kappa*Dt(i)*(X(i,:)-C(i,:))...
            + sqrt(2*D*Dt(i))*randn(1,2);
    end
end

%Final particle position
i=N;
if option.OU
    %use Ornstein-Uhlenbeck solution
    if z(i)==1
        X(i+1,:)= C(i,:)...
            + (X(i,:)-C(i,:))*exp(-kappa*Dt(i))...
            + sqrt((D/kappa)*(1-exp(-2*kappa*Dt(i))))*randn(1,2);
    else
        X(i+1,:)=X(i,:)+sqrt(2*D*Dt(i))*randn(1,2);
    end
else
    %Euler approximation to HPW model SDEs
    X(i+1,:)= X(i,:)...
        -(z(i))*kappa*Dt(i)*(X(i,:)-C(i,:))...
        + sqrt(2*D*Dt(i))*randn(1,2);
end
    
  

%give sample path as (xpos,ypos,time)
Y=[X [0; cumsum(Dt)]];

%add noise
if ~isempty(sigma)
    %store the true particle position
    X = Y;
    Traj.X = X;
    %simulate the observed particle position
    Y(:,1:2)= X(:,1:2) + sigma*randn(size(Y(:,1:2)));
    
    Traj.sigma = sigma;
end

Traj.Y=Y;
Traj.parameters=parameters;

Traj.z=z;
Traj.C=C;

Traj.type='Simulation';
Traj.Model='Harmonic Potential Well (HPW) Model';


end











