function [NbyTAbyCbyR, nTrials] = getNeuralDataTensor(Task, condIds)
%
% inputs:
%       Task: 
%       condIds: condition Ids from task table 
%
% outputs: 
%       NbyTAbyCbyR: [Neurons x Time x Conditions x Trials] data tensor.
%
% EMT 2019-12-30

N = Task.MU(condIds);

% build tensor over conditions
tensorTmp = zeros(0,0,0);  % initialize empty array of 1x3 to allow for first loop iteration
for iC = 1:length(condIds)
    tensorTmp = TensorUtils.catPad(4,tensorTmp, N(iC).rate);        
end

% strip off the zeros added for first loop iteration
tensorTmp(:,:,:,1) = []; 

% reshape for final output dimensions
NbyTAbyCbyR = permute(tensorTmp,[2 1 4 3]);

nTrials = [N(:).nTrials];