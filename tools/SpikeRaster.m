function [xPoints, yPoints] = SpikeRaster(spikes,varargin)
%   Inputs:
%       M x 1 cell of spike times:
%       M is the number of trials and each cell contains a 1 x N vector
%       of spike times. Units should be in seconds.
%
%   Output:
%       xPoints - vector of x points used for the plot.
%       yPoints - vector of y points used for the plot.
% 
%       MarkerFormat - default marker is a gray dot with size 1. Used for
%           scatter type plots only.
% 
%       AutoLabel - default 0.
%           Automatically labels x-axis as 'Time (ms)' or 'Time (s)' and
%           y-axis as 'Trial'.
% 
%       XLimForCell - default [NaN NaN].
%           Sets x-axis window limits if using cell spike time data. If
%           unchanged, the default limits will be 0.05% of the range. For
%           better performance, this parameter should be set.
% 
%       SpikeDuration - default 0.001 (1 millisecond).
%           Sets the horizontal spike length for cell spike time data.
%
%       RelSpikeStartTime - default 0 seconds.
%           Determines the starting point of the spike relative to the time
%           indicated by spike times or time bins. For example, a relative
%           spike start time of -0.0005 would center 1ms spikes for a
%           horzline plot of binary spike data.
%
%       rasterWindowOffset - default NaN
%           Exactly the same as relSpikeStartTime, but unlike
%           relSpikeStartTime, the name implies that it can be used to make
%           x-axis start at a certain time. If set, takes precedence over
%           relSpikeStartTime.
%  
%% Set Defaults and Load optional arguments
MarkerFormat.MarkerSize = 1;
MarkerFormat.Color = [0.2 0.2 0.2];
MarkerFormat.LineStyle = '.';

p = inputParser;
p.addRequired('spikes',@(x) iscell(x));
p.addParamValue('AutoLabel',1, @islogical);
p.addParamValue('XLimForCell',[NaN NaN],@(x) isnumeric(x) && isvector(x));
p.addParamValue('SpikeDuration',0.001,@(x) isnumeric(x) && isscalar(x));
p.addParamValue('RelSpikeStartTime',0,@(x) isnumeric(x) && isscalar(x));
p.addParamValue('RasterWindowOffset',NaN,@(x) isnumeric(x) && isscalar(x));
p.addParamValue('MarkerFormat',MarkerFormat,@isstruct);
p.parse(spikes,varargin{:});

spikes = p.Results.spikes;
autoLabel = p.Results.AutoLabel;
spikeDuration = p.Results.SpikeDuration;
xLimForCell = p.Results.XLimForCell;
relSpikeStartTime = p.Results.RelSpikeStartTime;
rasterWindowOffset = p.Results.RasterWindowOffset;

if ~isnan(rasterWindowOffset) && relSpikeStartTime==0
    relSpikeStartTime = rasterWindowOffset;
elseif ~isnan(rasterWindowOffset) && relSpikeStartTime~=0
    disp(['Warning: RasterWindoWOffset and RelSpikeStartTime perform the same function. '...
        'The value set in RasterWindowOffset will be used over RelSpikesStartTime']);
    relSpikeStartTime = rasterWindowOffset;
end
if ~isvector(spikes)
        error('Spike cell array must be an M x 1 vector.')
    end
    trialIsVector = cellfun(@isvector,spikes);
    if sum(trialIsVector) < length(spikes)
        error('Cells must contain 1 x N vectors of spike times.');
    end
    
    % Now make sure cell array is M x 1 and not 1 x M.
    if size(spikes,2) > 1 && size(spikes,1) == 1
        spikes = spikes';
    end
    
    % Make sure each trial is 1 x N and not N x 1
    nRowsInTrial = cellfun(@(x) size(x,1),spikes);
    % If there is more than 1 row in any trial, add a warning, and
    % transpose those trials. Allows for empty trials/cells (nRows > 1
    % instead of > 0).
    if sum(nRowsInTrial > 1) > 0 
        trialsToReformat = find(nRowsInTrial > 1);
        disp('Warning - some cells (trials) have more than 1 row. Those trials will be transposed.');
        for t = trialsToReformat
            spikes{trialsToReformat} = spikes{trialsToReformat}';
        end
    end
 %% 
    nTrials = length(spikes);
    
    % Find x-axis limits that aren't manually set (default [NaN NaN]), and
    % automatically set them. This is because we don't assume spikes start
    % at 0 - we can have negative spike times.
    limitsToSet = isnan(xLimForCell);
    if sum(limitsToSet) > 0
        % First find range of spike times
        minTimes = cellfun(@min,spikes,'UniformOutput',false);
        minTime = min( [ minTimes{:} ] );
        maxTimes = cellfun(@max,spikes,'UniformOutput',false);
        maxTime = max( [ maxTimes{:} ] );
        timeRange = maxTime - minTime;
        
        % Find 0.05% of the range.
        xStartOffset = relSpikeStartTime - 0.0005*timeRange;
        xEndOffset = relSpikeStartTime + 0.0005*timeRange + spikeDuration;
        newLim = [ minTime+xStartOffset, maxTime+xEndOffset ];
        xLimForCell(limitsToSet) = newLim(limitsToSet);
        % End result, if both limits are automatically set, is that the x
        % axis is expanded 0.1%, so you can see initial and final spikes.
    end
%     xlim(xLimForCell);
%     ylim([0 nTrials+1]);
    xPoints = [ spikes{:} ];
        
        % Getting the trials is trickier. 3 steps:
        % 1. First convert all the spike times into ones.
        trials = cellfun( @(x) {ones(size(x))}, spikes );
        % 2. Then multiply by trial number.
        for trialNum = 1:length(spikes)
            trials{trialNum} = trialNum*trials{trialNum};
        end
        % 3. Finally convert into a vector
        yPoints = [ trials{:} ];
        
end
