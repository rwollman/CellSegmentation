function [Ca,T,CaStk] = measureMCF10ACa2(MD,well,Lbl,varargin)

arg.timefunc = @(t) true(size(t)); 
arg = parseVarargin(varargin,arg); 


T = MD.getSpecificMetadata('TimestampFrame','Position',well,'Channel','Red','timefunc',arg.timefunc); 
T = cat(1,T{:}); 
CaStk = stkread(MD,'Position',well,'Channel','Red','timefunc',arg.timefunc);

%% register
CaStk = Lbl.Reg.register(CaStk,T); 

%% subtrack bacground (for all stack at once...); 
CaStk = backgroundSubtraction(CaStk); 
fprintf('Finsihed subtracting background')

%% do actual measurements
Ca = meanIntensityPerLabel(Lbl,CaStk,T,'func','median','type','nuc');

