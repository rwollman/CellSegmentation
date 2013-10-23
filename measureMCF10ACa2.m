function [Ca,T,CaStk] = measureMCF10ACa2(MD,well,Lbl,varargin)

T = MD.getSpecificMetadata('TimestampFrame','Position',well,'Channel','Red'); 
T = cat(1,T{:}); 
CaStk = stkread(MD,'Position',well,'Channel','Red'); 
% subtrack bacground (for all stack at once...); 
CaStk = backgroundSubtraction(CaStk,'timeskip',10); 
fprintf('Finsihed subtracting background')
Ca = meanIntensityPerLabel(Lbl,CaStk,T,'func','median','type','nuc');

