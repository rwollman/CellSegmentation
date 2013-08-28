function [Ca,T,CaStk] = measureMCF10ACa2(MD,well,Lbl,varargin)

T = MD.getSpecificMetadata('TimestampFrame','Position',well,'Channel','Red'); 
CaStk = stkread(MD,'Position',well,'Channel','Red'); 
% subtrack bacground (for all stack at once...); 
CaStk = backgroundSubtraction(CaStk); 
Ca = meanIntensityPerLabel(Lbl,CaStk,T,'func','median','type','nuc');

