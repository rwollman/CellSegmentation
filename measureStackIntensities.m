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
arg.samplingdensity = 15; 
arg.group_registration_by_acq = false; 
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




%% subtrack bacground (for all stack at once...); 
if arg.background
    switch arg.mskmethod
        case 'none'
            msk = true(size(measurementStack(:,:,1))); 
        case 'stack'
            ix = randi(numel(measurementStack),10000,1); 
            thrsh = prctile(measurementStack(ix),5); 
            mskStk = stkfun(@(m) m>thrsh,measurementStack);
            msk = mean(mskStk,3)==1; 
        case '5percentile'
            msk = nanmean(measurementStack,3);
            msk = msk>prctile(msk(:),5);
        case '5percentile_eroded'
            msk = nanmean(measurementStack,3);
            msk = msk>prctile(msk(:),5);
            msk = imerode(msk,strel('disk',50)); 
        case 'alpha-mean'
            prj = nanmean(measurementStack,3);
            p=prctile(prj(:),[5 95]);
            msk = false(size(prj)); 
            msk(prj>p(1) & prj<p(2))=true; 
            
        otherwise
            error('Mask call for background subtractin not supported!, check for typos...')
    end
    
    %% remove from the msk pixel that could belog to cells
    possiblecellbw = mean(Lbl.getLbls('base',Lbl.T),3)>0;
    msk(possiblecellbw)=false; 
    
    %% perform background subtraction
    measurementStack = backgroundSubtraction(measurementStack,'msk',msk,'smoothstk',false,'smoothmethod',arg.background_smooth, ...
                                                                         'percentile', arg.percentile,'samplingdensity',arg.samplingdensity);
    fprintf('Finsihed subtracting background\n')
end

%% register
if arg.register && isa(Lbl.Reg,'Registration') 
    if arg.group_registration_by_acq
        acq = MD.getSpecificMetadata('acq','Position',well,'Channel',arg.channel);
        [~,~,grouping] = unique(acq);
        measurementStack = Lbl.Reg.register(measurementStack,T,'grouping',grouping);
    else
        measurementStack = Lbl.Reg.register(measurementStack,T);
    end
    fprintf('Finished registration\n')
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
