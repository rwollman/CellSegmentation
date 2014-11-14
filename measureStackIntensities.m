function [Ca,T,CaStk] = measureStackIntensities(MD,well,Lbl,varargin)

arg.timefunc = @(t) true(size(t));
arg.channel = ''; 
arg.background = true; 
arg.register = true; 
arg.positiontype='Position';
arg.cellregiontouse = 'nuc'; 
arg.mskmethod = '5percentile';
arg.background_smooth = 'spline'; 

arg = parseVarargin(varargin,arg); 


if isempty(arg.channel)
    error('Channel is a requried argument!'); 
end


[CaStk,indx] = stkread(MD,arg.positiontype,well,'Channel',arg.channel,'timefunc',arg.timefunc,'sortby','TimestampFrame');
T = MD.getSpecificMetadataByIndex('TimestampFrame',indx); 

if iscell(T)
    CaStk = CaStk(:,:,~cellfun(@isempty,T));
    T = cat(1,T{~cellfun(@isempty,T)});
end


%% register
if arg.register && isa(Lbl.Reg,'Registration') 
    CaStk = Lbl.Reg.register(CaStk,T);
    fprintf('Finished registration')
end

%% subtrack bacground (for all stack at once...); 
if arg.background
    switch arg.mskmethod
        case 'none'
            msk = true(size(CaStk(:,:,1))); 
        case '5percentile'
            msk = nanmean(CaStk,3);
            msk = msk>prctile(msk(:),5);
        case '5percentile_eroded'
            msk = nanmean(CaStk,3);
            msk = msk>prctile(msk(:),5);
            msk = imerode(msk,strel('disk',50)); 
        otherwise
            error('Mask call for background subtractin not supported!, check for typos...')
    end
    CaStk = backgroundSubtraction(CaStk,'msk',msk,'smoothstk',false,'smoothmethod',arg.background_smooth);
    fprintf('Finsihed subtracting background')
end

%% do actual measurements
if iscell(arg.cellregiontouse)
    Ca = cell(size(arg.cellregiontouse));
    for i=1:numel(arg.cellregiontouse)
        Ca{i} = meanIntensityPerLabel(Lbl,CaStk,T,'func','mean','type',arg.cellregiontouse{i});
    end
else
    Ca = meanIntensityPerLabel(Lbl,CaStk,T,'func','mean','type',arg.cellregiontouse);
end
