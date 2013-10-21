function [Ca,Tca,CaStk] = measureMCF10ACa2(MD,well,Lbl,varargin)

%% input arguments
arg.cachannel = 'Red';
arg.positiontype = 'Position'; 
arg = parseVarargin(varargin,arg); 

%% read Ca stack
Tca = MD.getSpecificMetadata('TimestampFrame',arg.positiontype ,well,'Channel',arg.cachannel); 
CaStk = stkread(MD,arg.positiontype ,well,'Channel',arg.cachannel); 

Tca = cat(1,Tca{:});
[Tca ,ordr]= sort(Tca);
CaStk = CaStk(:,:,ordr); 

%% subtract  bacground (for all stack at once...); 
CaStk = backgroundSubtraction(CaStk);
% CaStk = backgroundSubtraction(CaStk,'timeskip',10); 


%% measure intensities
Ca = meanIntensityPerLabel(Lbl,CaStk,Tca,'func','median','type','nuc');

