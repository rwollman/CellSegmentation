function [Ca,Tca,CaStk] = measureMCF10ACa2(MD,well,Lbl,varargin)

%% input arguments
arg.cachannel = 'Red'; 
arg = parseVarargin(varargin,arg); 

%% read Ca stack
Tca = MD.getSpecificMetadata('TimestampFrame','Position',well,'Channel',arg.cachannel); 
CaStk = stkread(MD,'Position',well,'Channel','Red','TimestampFrame',Tca); 

Tca = cat(1,Tca{:});
[Tca ,ordr]= sort(Tca);
CaStk = CaStk(:,:,ordr); 

%% subtract  bacground (for all stack at once...); 
CaStk = backgroundSubtraction(CaStk); 

%% measure intensities
Ca = meanIntensityPerLabel(Lbl,CaStk,T,'func','median','type','nuc');

