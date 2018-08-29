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