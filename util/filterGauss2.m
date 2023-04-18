function OUT  = filterGauss2(IN, varargin)
% usage: out = filterGauss(in, SD);
%
% Filters an input data matrix with a gaussian smoothing kernel

% For convenience, the output is the same length as the input, but of 
% course the very early and late values are not 'trustworthy'.
%
% same filtering as in 'FilterSpikes' except the output is the same orientation as the input
% also more robust if SD is not a gaussian
% %
%
% EMT: updated on 2021-03-11 to generalize to arbitrary sampling rates


p = inputParser();
p.addRequired('IN',@isnumeric)  % data matrix
p.addParameter('gaussWidthMs',25, @isnumeric)  % SD in ms
p.addParameter('Fs',1000,@isnumeric) % Sampling frequency in seconds
% p.addParameter('gausswidth',8,@isnumeric) % number of SDs for extent of smoothin kernel

p.parse(IN, varargin{:})

gaussWidthMs = p.Results.gaussWidthMs;
Fs = p.Results.Fs;
% gausswidth = p.Results.gausswidth;

SDsamp = (Fs/1000)*gaussWidthMs;
SDrounded = 2 * round(SDsamp/2);  % often we just need SD to set the number of points we keep and so forth, so it needs to be an integer.
    

%% whatever orientation it came in, make it a row vector
% s = size(IN);
% 
% if s(1) > s(2)   % if it came in as a column vector
%     IN = IN';  % make it a row vector
%     flip = true;  % and remember to flip back at the very end
%     numV = s(2);
% else
%     flip = false;
%     numV = s(1);
% end

% convert to tall matrix
transposed = false;
if size(IN,2) > size(IN,1)
    IN = IN';
    transposed = true;
end


%% compute the normalized gaussian kernel
gaussWidthSamp = 8*SDsamp;
F = normpdf(1:gaussWidthSamp, gaussWidthSamp/2, SDsamp);
F = F/(sum(F));
    

outTmp = zeros(size(IN));
for ii = 1:size(IN,2)

    in = IN(:,ii);
    
    
    %% pad, filter
    shift = floor(length(F)/2); % this is the amount of time that must be added to the beginning and end;
    last = length(in); % length of incoming data
    prePad = makecol(zeros(1,shift)+mean(in(1:SDrounded)));
    postPad = makecol((zeros(1,shift)+mean(in(last-SDrounded:last))));
    preFilt = [prePad; in; postPad]; % pads the beginning and end with copies of first and last value (not zeros)
    postFilt = filter(F,1,preFilt); % filters the data with the impulse response in Filter
    
    %% trim
    
    postFiltTrim = postFilt(2*shift:length(postFilt)-1);  % Shifts the data back by 'shift', half the filter length
    outTmp(:,ii) = postFiltTrim;
       
end

% flip orientation if necessary, and send out
if transposed  % if it came in as a row vector
    OUT = outTmp';
else
    OUT = outTmp;
end


end


% function to ensure tall matrix regardless of input size
function vec = makecol( vec )

% transpose if it's currently a row vector (unless its 0 x 1, keep as is)
if (size(vec,2) > size(vec, 1) && isvector(vec)) && ~(size(vec, 1) == 0 && size(vec, 2) == 1)
    vec = vec';
end

if size(vec, 1) == 1 && size(vec, 2) == 0
    vec = vec';
end

end

