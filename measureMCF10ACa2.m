function [Ca,Tca,CaStk] = measureMCF10ACa2(MD,well,Lbl,varargin)


arg.cachannel='Red';
arg.positiontype='Position';
arg.timefunc = @(t) true(size(t)); 

arg = parseVarargin(varargin,arg); 


%% read Ca stack
T = MD.getSpecificMetadata('TimestampFrame',arg.positiontype,well,'Channel',arg.cachannel,'timefunc',arg.timefunc); 
T = cat(1,T{:}); 
CaStk = stkread(MD,arg.positiontype,well,'Channel',arg.cachannel,'timefunc',arg.timefunc);

%% subtract  bacground (for all stack at once...); 
CaStk = backgroundSubtraction(CaStk);



%% register
CaStk = Lbl.Reg.register(CaStk,T); 

%% subtrack bacground (for all stack at once...); 
CaStk = backgroundSubtraction(CaStk); 
fprintf('Finsihed subtracting background')

%% do actual measurements
Ca = meanIntensityPerLabel(Lbl,CaStk,T,'func','median','type','nuc');

