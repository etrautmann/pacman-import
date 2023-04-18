function [spikeTimesOut] = convertSpikeTimeIndices(spikeTimes, Fs1, Fs2)
% Convert a sparse matrix of spike times from one time base to another
% Anticipated use case is coverting a sparse matrix of spike times in units
% of samples on the IMEC headstage (sampled at 30kHz) into a sparse matrix
% of spike times in the time base of the NIDAQ (sampled at 32kHz). 
%
% inputs; spiketimes [nSample x nUnit] sparse array of spike times
%           Fs1 sampling frequency of origin time base (e.g. 30000 for IMEC headstage
%           Fs2 sampling frequency of output time base (e.g. 32000 for NIDAQ)
%
% %Note: for precision, both the IMEC headstage sampling clock rate and
% NIDAQ clock rate should be read in from the .meta files associated with
% each, not assumed to be perfectly accurate
% EMT 2021-03-24


[inds, cols] = find(spikeTimes);

inds2 = round(inds*Fs2/Fs1);

spikeTimesOut = sparse(inds2, cols, true);