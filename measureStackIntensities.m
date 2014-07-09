function [Ca,T,CaStk] = measureStackIntensities(MD,well,Lbl,varargin)

arg.timefunc = @(t) true(size(t));
arg.channel = ''; 
arg.background = true; 
arg.register = true; 
arg.region = 'nuc'; 
arg = parseVarargin(varargin,arg); 

if isempty(arg.channel)
    error('Channel is a requried argument!'); 
end


[CaStk,indx] = stkread(MD,'Position',well,'Channel',arg.channel,'timefunc',arg.timefunc,'sortby','TimestampFrame');
T = MD.getSpecificMetadataByIndex('TimestampFrame',indx); 
T = cat(1,T{:}); 

%% register
if arg.register
    CaStk = Lbl.Reg.register(CaStk,T);
    fprintf('Finished registration')
end

%% subtrack bacground (for all stack at once...); 
if arg.background
    msk = nanmean(CaStk,3);
    msk = msk>prctile(msk(:),5);
    CaStk = backgroundSubtraction(CaStk,'msk',msk);
    fprintf('Finsihed subtracting background')
end

%% do actual measurements
Ca = meanIntensityPerLabel(Lbl,CaStk,T,'func','mean','type',arg.region);

