function [PerCellRatio,T,RatioStk] = measureStackRatio(MD,Lbl,well,varargin)

arg.channels ={}; 
arg.interpolatetochanel = ''; 

arg.numweights
arg.denomweights

arg.numeratorordenominator = true(0,1); % true / false if numerator (true) or denominator (false)
arg.weights

arg.numeratorstack = {}; 
arg.numeratorstackweights = []; 
arg.denominatorstacks = {}; 
arg.denominatorstackweights = []; 


%% read stacks
Channels = arg.channels; 
Stacks = cell(size(arg.channels)); 
Timevectors = cell(size(Stacks)); 

for i=1:numel(Stacks)
    %% read stack
    
    %% register stack
    
end

%% interpolate if needed
if ~isempty(arg.interpolatetochanel)
    [~,ix_channelstointerpolate] = setdiff(Channels,arg.interpolatetochanel)
    for i=1:numel(ix_channelstointerpolate)
    %% interploate stack
    
    end
end

sz = size(Stacks{1}); % they should all have the same size

%% backgroud subtraction 

%% create ratio stack
Numer = zeros(sz); 
Denomi = zeros(sz); 

for i=1:numel(Channels)
    Numer = Numer + arg.numweights(i)*Stacks{i}; 
    Denomi = Denomi + arg.denomweights(i)*Stacks{i}; 
end

ratio = Numer./ Denomi; 

%%  measure per label 
PerCellRatio = meanIntensityPerLabel(Lbl,ratio,T,'func','median','type','cyto'); 

%% make Ratio stack to allow for nice visualization
if nargout ==3
    f=imresize(ratio(:,:,1),0.333);
    RatioStk = zeros(size(f,1),size(f,2),size(ratio,3),'single');
    Imin = prctile(ratio(unidrnd(numel(ratio),10000,1)),5); 
    for i=1:size(ratio,3)
        lbl =getLbls(Lbl,'cyto',T(i));
        lbl = imresize(lbl,0.333);
        f=imresize(ratio(:,:,i),0.333);
        f(~lbl)=Imin;
        RatioStk(:,:,i)=f;
    end
end