function [measurement,T,measurementStack, cellids] = measureStackIntensities(MD,well,Lbl,varargin)

arg.timefunc = @(t) true(size(t));
arg.channel = ''; 
arg.background = true; 
arg.register = true; 
arg.positiontype='Position';
arg.cellregiontouse = 'nuc'; 
arg.mskmethod = '5percentile';
arg.background_smooth = 'spline';
arg.func = 'mean';
arg.percentile = 0; % Parameter passed to backgroundSubtraction

arg = parseVarargin(varargin,arg); 


if isempty(arg.channel)
    error('Channel is a requried argument!'); 
end


[measurementStack,indx] = stkread(MD,arg.positiontype,well,'Channel',arg.channel,'timefunc',arg.timefunc,'sortby','TimestampFrame');
T = MD.getSpecificMetadataByIndex('TimestampFrame',indx); 

if iscell(T)
    measurementStack = measurementStack(:,:,~cellfun(@isempty,T));
    T = cat(1,T{~cellfun(@isempty,T)});
end


%% register
if arg.register && isa(Lbl.Reg,'Registration') 
    measurementStack = Lbl.Reg.register(measurementStack,T);
    fprintf('Finished registration\n')
end

%% subtrack bacground (for all stack at once...); 
if arg.background
    switch arg.mskmethod
        case 'none'
            msk = true(size(measurementStack(:,:,1))); 
        case '5percentile'
            msk = nanmean(measurementStack,3);
            msk = msk>prctile(msk(:),5);
        case '5percentile_eroded'
            msk = nanmean(measurementStack,3);
            msk = msk>prctile(msk(:),5);
            msk = imerode(msk,strel('disk',50)); 
        otherwise
            error('Mask call for background subtractin not supported!, check for typos...')
    end
    measurementStack = backgroundSubtraction(measurementStack,'msk',msk,'smoothstk',false,'smoothmethod',arg.background_smooth, 'percentile', arg.percentile);
    fprintf('Finsihed subtracting background\n')
end

%% do actual measurements
cellids={}; 
if iscell(arg.cellregiontouse)
    measurement = cell(size(arg.cellregiontouse));
    cellids = cell(size(arg.cellregiontouse));
    for i=1:numel(arg.cellregiontouse)
        measurement{i} = meanIntensityPerLabel(Lbl,measurementStack,T,'func',arg.func,'type',arg.cellregiontouse{i});
%         cellids{i} = Lbl.cellids;
    end
else
    measurement = meanIntensityPerLabel(Lbl,measurementStack,T,'func',arg.func,'type',arg.cellregiontouse);
%     cellids = Lbl.cellids;
end
