function [zProp,LogzDensity,LogzPropDensity,PropType] = MetropolisHastingsForHiddenState(z,options)
%Metropolis-Hastings proposal for a binary hidden (optionally Markov) state model

%Given a current hidden state sequence, z, calculates a proposed hidden state
%sequence using the proposal density q.
%Also calculates the densities q(z->z'), q(z'->z).
%So that the acceptance, \alpha(z->z')=(p(z')*q(z'->z))/(p(z)*q(z'->z))
%can be calculated, where p(z) is the density proportional to the density
%of interest.

%There are a number of options for q
%Optional proposal moves:
%A. Blocked move
%B. Shift move

%For blocked move, there are three options
%1. Move entire block to state 0
%2. Move entire block to state 1
%3. Simulate a Markov chain for the new block states (necessary for
%reversibility of the proposal distribution)

% INPUT:
% z - current hidden state sequence
% options.ProposalOptions - choose how to propose a new hidden state
% sequence, choices are: "Blocked" and "Shift"
% options.StateNames - the names of the hidden states
% options.BlockToMove - optional - specify the block to move rather than
% sampling it
% options.MaxShift - maximum shift size for the shift move

% OUTPUT:
% zProp - the proposed hidden state sequence
% LogzDensity - proposal probablity q(z'->z), i.e. probability of proposing
% a move from z' to z under the proposal choice
% LogzPropDensity - proposal probablity q(z->z')
% PropType - the proposal choice
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


N=length(z);

%rename states as 0,1 for simplicity 
%check if there are only 0,1 in z
if ~((sum(z==0)+sum(z==1))==N)
    z(z==options.StateNames(1))=0;
    z(z==options.StateNames(2))=1;
end

%options are any combination of "Blocked" and "Shift"
ProposalOptions=options.ProposalOptions;
%if both "Blocked" and "Shift" are selected the randomly (uniform dist)
%select between them
if iscell(ProposalOptions)
    PropType=ProposalOptions{randi(length(ProposalOptions))};
else
    PropType=options.ProposalOptions;
end

switch PropType    
    case 'Blocked'
        %
        %Blocked move
        %
        
        %transition probabilities (can these can be tuned during burn in?)
        p_01=0.01;
        p_10=0.01;
                                
        zProp=z;
        
        %update a single block
        %sample block size
        if isfield(options,'BlockToMove')
            BlockToMove=options.BlockToMove;
        else
            BlockToMove = SampleBlockToMove(N,options);                        
        end
        
        BlockLength=length(BlockToMove);
        
        
        %propose z by moving a block, length=BlockLength
        %either all zProp=1,all zProp=0, or simulate Markov chain
        %also calculate proposal density for z', q(z')
        UnifRand=rand;
        if UnifRand <1/3
            zProp(BlockToMove)=1;
            
            ...also possible to simulate all ones, so need extra term in proposal density
                i=BlockToMove(1);
            if i==1
                LogzPropDensity=log(1/3 ...
                    +(1/3)*zProp(i)*(p_01/(p_10+p_01))...
                    *(1-p_10)^(BlockLength-1));
            else
                LogzPropDensity=log(1/3 ...
                    +(1/3)*(...
                    zProp(i)*zProp(i-1)*(1-p_10)...
                    +zProp(i)*(1-zProp(i-1))*p_01...
                    )...
                    *(1-p_10)^(BlockLength-1));
            end
        elseif UnifRand < 2/3
            zProp(BlockToMove)=0;
            ...also possible to simulate all zeros, so need extra term in proposal density
                i=BlockToMove(1);
            if i==1
                LogzPropDensity=log(1/3 ...
                    +(1/3)*(1-zProp(i))*(p_10/(p_10+p_01))...
                    *(1-p_01)^(BlockLength-1));
            else
                LogzPropDensity=log(1/3 ...
                    +(1/3)*(...
                    (1-zProp(i))*zProp(i-1)*p_10...
                    +(1-zProp(i))*(1-zProp(i-1))*(1-p_01)...
                    )...
                    *(1-p_01)^(BlockLength-1));
            end
        else
            %simulate Markov chain
            zProp(BlockToMove)=0;
            zProposalDensity=zeros(BlockLength,1);
            j=1;
            for i=BlockToMove
                if i==1
                    %if first timepoint simulate from stationary
                    %distributions
                    zProp(i)=rand<p_01/(p_10+p_01);
                    
                    zProposalDensity(j)=...
                        zProp(i)*(p_01/(p_10+p_01))...
                        +(1-zProp(i))*(p_10/(p_10+p_01));
                    
                    j=j+1;
                else
                    zProp(i)=rand<(zProp(i-1)*(1-p_10)+(1-zProp(i-1))*p_01);
                    
                    zProposalDensity(j)=...
                        zProp(i)*zProp(i-1)*(1-p_10)...
                        +zProp(i)*(1-zProp(i-1))*p_01...
                        +(1-zProp(i))*zProp(i-1)*p_10...
                        +(1-zProp(i))*(1-zProp(i-1))*(1-p_01);
                    
                    j=j+1;
                end
            end
            
            if all(zProp(BlockToMove)==1) || all(zProp(BlockToMove)==0)
                zProposalDensity=1/3+1/3*prod(zProposalDensity);
            else
                zProposalDensity=1/3*prod(zProposalDensity);
            end
            LogzPropDensity=log(zProposalDensity);
        end
        
        %calculate proposal density for current z, q(z) given zProp
        if all(z(BlockToMove)==1)
            i=BlockToMove(1);
            if i==1
                LogzDensity=log(1/3 ...
                    +(1/3)*z(i)*(p_01/(p_10+p_01))...
                    *(1-p_10)^(BlockLength-1));
            else
                LogzDensity=log(1/3 ...
                    +(1/3)*(...
                    z(i)*zProp(i-1)*(1-p_10)...
                    +z(i)*(1-zProp(i-1))*p_01...
                    )...
                    *(1-p_10)^(BlockLength-1));
                
            end
        elseif all(z(BlockToMove)==0)
            i=BlockToMove(1);
            if i==1
                LogzDensity=log(1/3 ...
                    +(1/3)*(1-z(i))*(p_10/(p_10+p_01))...
                    *(1-p_01)^(BlockLength-1));
            else
                LogzDensity=log(1/3 ...
                    +(1/3)*(...
                    (1-z(i))*zProp(i-1)*(p_10)...
                    +(1-z(i))*(1-zProp(i-1))*(1-p_01)...
                    )...
                    *(1-p_01)^(BlockLength-1));
            end
        else
            zDensity=zeros(BlockLength,1);
            j=1;
            for i=BlockToMove
                if i==1
                    zDensity(j)=...
                        z(i)*(p_01/(p_10+p_01))...
                        +(1-z(i))*(p_10/(p_10+p_01));
                    
                    j=j+1;
                else
                    zDensity(j)=...
                        z(i)*z(i-1)*(1-p_10)...
                        +z(i)*(1-z(i-1))*p_01...
                        +(1-z(i))*z(i-1)*p_10...
                        +(1-z(i))*(1-z(i-1))*(1-p_01);
                    
                    j=j+1;
                end
            end
            zDensity=(1/3)*prod(zDensity);
            LogzDensity=log(zDensity);
        end
        
        
        
        
        
    case 'Shift'
        %Shift move
        
        %pick switch event to shift
        SwitchEvents=find(z(1:end-1)~=z(2:end));
        zProp=z;
            
        if ~isempty(SwitchEvents)
            if length(SwitchEvents)==1
                k=SwitchEvents;
            else
                k=randsample(SwitchEvents,1);
                while k==1 || k==N
                    k=randsample(SwitchEvents,1);
                end
            end
            MaxShift=options.MaxShift;
            %find nearest switch to the switch to move
            NearestSwitch=min(k-SwitchEvents(SwitchEvents~=k));
            N_max=floor(0.5*(min([abs(k-NearestSwitch) ...
                k-1 ...
                N-k...
                MaxShift])));                    
            
            if N_max>0
                n=randi(N_max);
                
                if rand<1/2
                    zProp(k:k+n)=zProp(k-1);
                    PropShift=k+n;
                else
                    zProp(k-n:k)=zProp(k+1);
                    PropShift=k-n;
                end
                
                
                %also calculate proposal density for z', q(z'|z)
                LogzDensity=1/(2*N_max);
                
                
                %calculate proposal density for current z, q(z|z') (given
                %zProp)
                PropSwitchEvents=find(zProp(1:end-1)~=zProp(2:end));
                
                %NearestSwitch is also the nearest switch to PropShift
                
                PropN_max=floor(0.5*(min([abs(PropShift-NearestSwitch) ...
                    PropShift-1 ...
                    N-PropShift...
                    MaxShift])));
                
                LogzPropDensity=1/(2*PropN_max);
                
            else
                LogzDensity=0;
                LogzPropDensity=0;
            end
        else %if no switch events
            LogzDensity=0;
            LogzPropDensity=0;
        end
end




if ~((sum(z==0)+sum(z==1))==N)
    %rename states as 0,1 for simplicity
    zPropTemp=zeros(N,1);
    zPropTemp(zProp==0)=options.StateNames(1);
    zPropTemp(zProp==1)=options.StateNames(2);
    
    zProp=zPropTemp;
end






end