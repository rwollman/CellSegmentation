function [Erk,T,FretStk] = measureEKARevFRET(MD,well,Lbl,varargin)
% Performs all the necessary image processing and calculate FRET ratio
% The code does the following processing before measurement:  
% 1. Image registration - using normxcorr2 on central croped image. works pretty well. 
% 2. Background subtraction - using a exponential fit for a min of central
%    image. 
% 3. bleedthrough correction for CFP and YFP (even if YFP is more sparsly
%    acquired) based on bleedthrough calibration. 


%% input arguments
arg.timefunc = @(t) true(size(t)); 
arg.cfp2fretbleedthrough = 0.95; 
arg.yfp2fretbleedthrough = 0.0284; 
arg.backgroundfilter = sum(fspecial('gauss',5,3));

arg.positiontype = 'Position'; 

arg.crop = [680 680 680 680]; 
arg.positiontype = 'Position';

arg.register = []; 

arg = parseVarargin(varargin,arg); 


%% get timepoints for measurements. 
T = MD.getSpecificMetadata('TimestampFrame','Channel','Cyan',arg.positiontype ,well,'timefunc',arg.timefunc);
T=cat(1,T{:}); 

%% Read Cyan / CyanToYellow stack and process them. 
c2y = stkread(MD,arg.positiontype ,well,'Channel','CyanToYellow','timefunc',arg.timefunc);
cfp = stkread(MD,arg.positiontype ,well,'Channel','Cyan','timefunc',arg.timefunc);

assert(~isempty(cfp),'Error CFP is empty - check argument calls'); 

% make sure that numel(T) is the same as (it should be but if not, try to
% correct as best as possible. 
if size(c2y,3)~=numel(T)
    T=linspace(min(T),max(T),size(c2y,3));
    T=T(:); 
end

%% register the stacks 
if isa(arg.register,'Registration')
% If arg.register exist than just use it
    Reg = arg.register;
    cfp = Reg.register(cfp,T); 
    c2y = Reg.register(c2y,T); 
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
    ref2d = imref2d([size(cfp,1) size(cfp,2)]);
    
    % repeat registration
    for i=1:numel(unq)
        ix=find(refix==unq(i));
        [~,mi] = min(minD(ix));
        [c2y(:,:,ix),Tfrms] = registerStack(c2y(:,:,ix),'reference',mi,'crop',arg.crop,'method','xcorr');
        for j=1:numel(ix)
            cfp(:,:,ix(j))=imwarp(cfp(:,:,ix(j)),Tfrms(j),'OutputView',ref2d);
        end
    end
end

%% measure YFP and interpolate for missing values
Tlbl = Lbl.T; 
Tlbl = Tlbl(arg.timefunc(Tlbl)); 
yfpsml = stkread(MD,arg.positiontype ,well,'Channel','Yellow','timefunc',arg.timefunc,'TimestampFrame',Tlbl);
assert(numel(Tlbl)==size(yfpsml,3),'The number of YFP slices is different then the timepoints in Lbl object'); 


% for yfp just subtract background without filtering. 
for i=1:size(yfpsml,3)
    m=imcrop(yfpsml(:,:,i),arg.crop); 
    yfpsml(:,:,i)=yfpsml(:,:,i)-min(m(:)); 
end

if size(yfpsml,3)==1
    yfp = repmat(yfpsml,[1 1 numel(T)]); 
elseif size(yfpsml,3)==size(c2y,3)
    yfp=yfpsml;
else
    yfprow = reshape(yfpsml,size(cfp,1)*size(cfp,2),size(yfpsml,3))';
    yfp = interp1(Tlbl,yfprow,T,'nearest','extrap');
    yfp = reshape(yfp',[Lbl.sz numel(T)]);
end
if isa(arg.register,'Registration')
    yfp = Reg.register(yfp,T); 
end
% save memory
clear yfpsml; 


%% perform background subtraction
cfp_bck=nan(size(cfp,3),1); 
c2y_bck=nan(size(c2y,3),1); 
for i=1:numel(c2y_bck)
    m1=imcrop(cfp(:,:,i),arg.crop); 
    cfp_bck(i)=min(m1(:)); 
    m2=imcrop(c2y(:,:,i),arg.crop); 
    c2y_bck(i)=min(m2(:)); 
end
% 

Tsec = (T(:)-min(T))*24*60;
n=numel(Tsec);
b=regress(log(c2y_bck(1:n)),[ones(n,1) Tsec(1:n)]); 
c2y_fit = exp(b(1))*exp(b(2)*Tsec); 

b=regress(log(cfp_bck(1:n)),[ones(n,1) Tsec(1:n)]); 
cfp_fit = exp(b(1))*exp(b(2)*Tsec);

for i=1:n
    c2y(:,:,i) = c2y(:,:,i)-c2y_fit(i); 
    cfp(:,:,i) = cfp(:,:,i)-cfp_fit(i); 
end

cfp=cfp(:,:,1:n); 
c2y = c2y(:,:,1:n); 

%% Bleedthrough calculations
ratio = (c2y-arg.cfp2fretbleedthrough*cfp-arg.yfp2fretbleedthrough*yfp)./cfp; 
clear c2y yfp cfp % save memory
% next two lines should do nothing, in some cases due to messed acqusition
% they are needed
[T,ordr]=sort(T); 
ratio = ratio(:,:,ordr); 

%% actual measurement 
Erk = meanIntensityPerLabel(Lbl,ratio,T,'func','median','type','cyto'); 

%% make FretStak to allow for nice visualization
if nargout ==3
    f=imresize(ratio(:,:,1),0.333);
    FretStk = zeros(size(f,1),size(f,2),size(ratio,3),'single');
    Imin = prctile(ratio(unidrnd(numel(ratio),10000,1)),5); 
    for i=1:size(ratio,3)
        lbl =getLbls(Lbl,'cyto',T(i));
        lbl = imresize(lbl,0.333);
        f=imresize(ratio(:,:,i),0.333);
        f(~lbl)=Imin;
        FretStk(:,:,i)=f;
    end
end

