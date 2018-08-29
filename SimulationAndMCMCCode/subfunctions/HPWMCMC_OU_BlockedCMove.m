function [C, LogCDensity]=HPWMCMC_OU_BlockedCMove(Traj,C,z,parameters,prior,options)
%if options.sample is 0 then don't actually do the move - just work out the
%Gibbs density

D=parameters(1);
D_C=parameters(2);
kappa=parameters(3);
D_est=parameters(4);


X=Traj.Y(:,1:2);
dX=diff(X);
Dt=diff(Traj.Y(:,3));
N=length(Dt);

%update a single block
%sample block size
if isfield(options,'BlockToMove')
    BlockToMove=options.BlockToMove;
else
    BlockToMove = SampleBlockToMove(N,options);
end
BlockSize=length(BlockToMove);

sigma_diag=zeros(N,1);
b=zeros(N,2);

for i=BlockToMove   
    tau1=(z(i)*(1-exp(-kappa*Dt(i)))^2)/...
        (D*(1-exp(-2*kappa*Dt(i)))/kappa);
    
    mu1=z(i)*(dX(i,:)+X(i,:)*(1-exp(-kappa*Dt(i)))/...
        (1-exp(-kappa*Dt(i))));
    
    if i~=1
        tau2=1/(2*Dt(i-1)*(D_C*z(i-1)+D_est*(1-z(i-1))));
        mu2=C(i-1,:);
    else       
        tau2=1/prior.sigma_C;
        mu2=prior.mu_C;
    end
    if i~=N
        tau3=1/(2*Dt(i)*(D_C*z(i)+D_est*(1-z(i))));
        mu3=C(i+1,:);
    else
        tau3=0;
        mu3=0;
    end
                               
    %diagonal values of covariance matrix sigma
    sigma_diag(i)=tau1+tau2+tau3;
         
    %CPrecisionMatrix*mu=b
    if i==BlockToMove(1)
        b(i,:)=tau1*mu1 + mu2*tau2;
    elseif i==BlockToMove(end)
        b(i,:)=tau1*mu1 + mu3*tau3;
    else
        b(i,:)=tau1*mu1;
    end

end
    
sigma_diag=sigma_diag(BlockToMove);
b=b(BlockToMove,:);

%Precision matrix blocked move
CPrecisionMatrix = sparse(BlockSize,BlockSize);

%diagonal elements           
CPrecisionMatrix(1:(BlockSize+1):end) = sigma_diag;
%C_{i,i+1} coefficients
CPrecisionMatrix((BlockSize+1):(BlockSize+1):end) = -1 ...
    ./(2.*Dt(BlockToMove(1:end-1)).*(D_est.*(1-z(BlockToMove(1:end-1)))+D_C.*z(BlockToMove(1:end-1))));
%C_{i,i-1} coefficients
CPrecisionMatrix(2:(BlockSize+1):end) = -1 ...
    ./(2.*Dt(BlockToMove(1:end-1)).*(D_est.*(1-z(BlockToMove(1:end-1)))+D_C.*z(BlockToMove(1:end-1))));

%Cholsky decomposition of precision matrix
[CPrecisionChol,~] = chol(CPrecisionMatrix);


%calculate mean vector mu_C
mu_C=CPrecisionMatrix\b;

if options.sample
    %sample mvn (without having to invert precision matrix)
    C(BlockToMove,:)=CPrecisionChol\randn(BlockSize,2) + mu_C;
    
    %calculate the MVN density for the sampled C
    LogCDensity=zeros(2,1);
    for i=1:2
        LogCDensity(i)=-(BlockSize/2)*log(2*pi)...
            +0.5*2*sum(log(diag(chol(CPrecisionMatrix))))...
            -0.5*(C(BlockToMove,i)-mu_C(:,i))'...
            *CPrecisionMatrix...
            *(C(BlockToMove,i)-mu_C(:,i));
    end
    LogCDensity=sum(LogCDensity);
else
    %just calculate the MVN density for C provided in options
    %(i.e. the probability of drawing options.CForDensity, given function input C, z
    
    LogCDensity=zeros(2,1);
    for i=1:2
        LogCDensity(i)=-(BlockSize/2)*log(2*pi)...
            +0.5*2*sum(log(diag(chol(CPrecisionMatrix))))...
            -0.5*(options.CForDensity(BlockToMove,i)-mu_C(:,i))'...
            *CPrecisionMatrix...
            *(options.CForDensity(BlockToMove,i)-mu_C(:,i));
    end
    LogCDensity=sum(LogCDensity);
end


end