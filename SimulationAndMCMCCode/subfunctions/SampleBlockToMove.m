function Blocks = SampleBlockToMove(N,options)
% Randomly split sequence of length N into blocks

% INPUT: N - length of the hidden state sequence
% options.MinBlockSize, options.MaxBlockSize - min and max sizes of the
% blocks
% options.MultipleBlocks - if=1 returns start and endpoints of multiple blocks
% (i.e. splits the whole sequence into blocks)
% if=0 returns all the timepoints for a single block to move
%
% OUTPUT: Blocks - the blocks to attempt MCMC moves on
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



%sample block size (same even for multiple blocks)
BlockSize=randi([options.MinBlockSize options.MaxBlockSize]);

if isfield(options,'MultipleBlocks') && options.MultipleBlocks
    BlockStarts=(1:BlockSize:(N-1));
    BlockEnds=BlockStarts+BlockSize;
    BlockEnds(BlockEnds>N)=N;
    Blocks=[BlockStarts' BlockEnds'];
else
    %sample the first timepoint in the block
    BlockStart=randi([1 N-BlockSize]);
    
    Blocks=BlockStart:BlockStart+BlockSize-1;
end


end