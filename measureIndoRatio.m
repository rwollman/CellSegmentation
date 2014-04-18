function [Ca,T,IndoStk] = measureIndoRatio(MD,well,Lbl,varargin)
% Performs all the necessary image processing and calculate FRET ratio
% The code does the following processing before measurement:  
% 1. Image registration - using normxcorr2 on central croped image. works pretty well. 
% 2. Background subtraction - using a exponential fit for a min of central
%    image. 
% 3. bleedthrough correction for CFP and YFP (even if YFP is more sparsly
%    acquired) based on bleedthrough calibration. 


%% input arguments
arg.timefunc = @(t) true(size(t)); 

arg.backgroundfilter = sum(fspecial('gauss',5,3));

arg.positiontype = 'Position'; 

arg.crop = [680 680 680 680]; 

arg.register = []; 
arg.background = true; 
arg = parseVarargin(varargin,arg); 


%% get timepoints for measurements. 
T = MD.getSpecificMetadata('TimestampFrame','Channel','IndoCaFree',arg.positiontype ,well,'timefunc',arg.timefunc);
T=cat(1,T{:}); 

%% Read CaBound / CaFree stack and process them. 
CaBound = stkread(MD,arg.positiontype ,well,'Channel','IndoCaBound','timefunc',arg.timefunc);
CaFree= stkread(MD,arg.positiontype ,well,'Channel','IndoCaFree','timefunc',arg.timefunc);

% make sure that numel(T) is the same as (it should be but if not, try to
% correct as best as possible. 
if size(CaBound,3)~=numel(T)
    T=linspace(min(T),max(T),size(CaBound,3));
    T=T(:); 
end

%% register the stacks 
if isa(arg.register,'Registration')
% If arg.register exist than just use it
    Reg = arg.register;
    CaFree = Reg.register(CaFree,T); 
    CaBound = Reg.register(CaBound,T); 
else
% registration happens such that every image will be registered to the
% reference image that was taken at the same time as the image used for
% labeling. This means that we will run the registerStack functions as many
% times as we have labels.

    % figure out who "belongs" to whon
    Tlbl = Lbl.T;
    Tlbl(~arg.timefunc(Tlbl))=[]; % remove Tlbl according to timefunc
    D = abs(repmat(Tlbl(:)',size(T,1),1)-repmat(T(:),1,numel(Tlbl)));
    [minD,refix]=min(D,[],2);
    unq = unique(refix);
    ref2d = imref2d([size(CaFree,1) size(CaFree,2)]);
    
    % repeat registration
    for i=1:numel(unq)
        ix=find(refix==unq(i));
        [~,mi] = min(minD(ix));
        [CaBound(:,:,ix),Tfrms] = registerStack(CaBound(:,:,ix),'reference',mi,'crop',arg.crop,'method','xcorr');
        for j=1:numel(ix)
            CaFree(:,:,ix(j))=imwarp(CaFree(:,:,ix(j)),Tfrms(j),'OutputView',ref2d);
        end
    end
end




%% perform background subtraction
if arg.background
    msk = nanmean(CaBound,3);
    msk = msk>prctile(msk(:),5);
    CaBound = backgroundSubtraction(CaBound,'msk',msk,'smoothstk',false);
    CaFree = backgroundSubtraction(CaFree,'msk',msk,'smoothstk',false);
end



%% Calculate ratio calculations
ratio = (CaBound./CaFree); 
clear CaBound CaFree % save memory

%% actual measurement 

Ca = meanIntensityPerLabel(Lbl,ratio,T,'func','median','type','nuc'); 

%% make IndoStak to allow for nice visualization
if nargout ==3
    f=imresize(ratio(:,:,1),0.333);
    IndoStk = zeros(size(f,1),size(f,2),size(ratio,3),'single');
    Imin = prctile(ratio(unidrnd(numel(ratio),10000,1)),5); 
    for i=1:size(ratio,3)
        lbl =getLbls(Lbl,'base',T(i));
        lbl = imresize(lbl,0.333);
        f=imresize(ratio(:,:,i),0.333);
        f(~lbl)=Imin;
        IndoStk(:,:,i)=f;
    end
end

