function [Lbl,NucLabels,CellLabels,CytoLabels,Timing] = segmentNucleiAndRingsWithoutTracking(MD,well,varargin)

arg.nuc_channel = 'DeepBlue'; 
arg.nuc_smooth1 = 5; % sigma of filtering done to improve thresholding. The size of the nuclei will be eroded by this sigma as well
arg.nuc_smooth2 = fspecial('gauss',15,9);
arg.nuc_smooth = strel('disk',5); 
arg.cyto_channels = {'Cyan','Yellow','Red'}; 
arg.positiontype = 'Position'; 
arg.register = []; % optional registration object
arg.timefunc = @(t) true(size(t));
arg.cyto_ringstrel = strel('disk',15); 
arg.cyto_transform='log'; % other alternative: linear
arg.cyto_thresholdmethod = 'otsu'; % other alternatives: gm,minerr,robust,localotsu,kmeans
arg.sz =  [2048  2064]; 
arg = parseVarargin(varargin,arg); 
    
t0=now; 

%% get timepoints for the Label matrices
T = MD.getSpecificMetadata('TimestampFrame','Channel',arg.nuc_channel,arg.positiontype,well,'timefunc',arg.timefunc);
T = cat(1,T{:});
nuc = stkread(MD,'TimestampFrame',T,'Channel',arg.nuc_channel); 
[T,ordr]= sort(T);
nuc=nuc(:,:,ordr); 

Timing.readnucimages=now-t0;
t0=now; 

%% register if requested or if Registration object was passed as an input arguemtn
if islogical(arg.register) && arg.register
    [nuc,Tforms] = registerStack(nuc,'crop',[400 400 1200 1200],'method','xcorr','reference',ceil(size(nuc,3)/3),'maxdisp',100);
    Reg = Registration(T,Tforms);
    arg.register = Reg; 
elseif ~isempty(arg.register) && isa(arg.register,'Registration')
    Reg = arg.register; 
    nuc = Reg.register(nuc,T); 
else
    Reg = []; 
end

Timing.registration=now-t0;
t0=now; 

%% first subtrack background for entire stack
% subtract bacgkround using a mask to avoid corner issues
nucprj = mean(nuc,3);
msk = nucprj>prctile(nucprj(:),5);
msk = imfill(msk,'holes');
msk  = bwareaopen(msk,10000);

nuc = backgroundSubtraction(nuc,'msk',logical(msk),'smoothstk',false);
NucLabels = zeros([arg.sz size(nuc,3)],'uint16');
CellLabels = zeros([arg.sz size(nuc,3)],'uint16');
CytoLabels = zeros([arg.sz size(nuc,3)],'uint16');

Timing.backgroundsubtraction=now-t0;
t0=now; 


%% for each well, segment all nuclei
nuc_smooth1 = arg.nuc_smooth1; 
nuc_smooth2=arg.nuc_smooth2; 
parfor i=1:numel(T)
    %% segment all nuclei
    nucfilt = imfilter(nuc(:,:,i),fspecial('gauss',2*nuc_smooth1,nuc_smooth1)) ;
    bw = optThreshold(nucfilt,'msk',msk,'method','localotsu','transform','none');
    pks = imregionalmax(imfilter(nucfilt,nuc_smooth2));
    bw  = bw &~bwareaopen(bw,1000);
    NucLabels(:,:,i) = segmentUsingSeeds(bw,bwlabel(pks & bw));
end

Timing.nucleisegmentation=now-t0;
t0=now; 

%% create a whole cell binary mask
BW = false(size(nuc)); 
for j=1:numel(arg.cyto_channels)
    MAPK = stkread(MD,'TimestampFrame',T,'Channel',arg.cyto_channels{j});
    if islogical(arg.register) && arg.register
        Reg.register(MAPK,T); 
    end
    MAPK = imfilter(MAPK,fspecial('gauss',5,3)); 
    parfor i=1:numel(T)
        %% segment cytoplasm
        BW(:,:,i) = BW(:,:,i) | optThreshold(MAPK(:,:,i),'msk',msk,'method',arg.cyto_thresholdmethod,'transform',arg.cyto_transform);
    end
end

Timing.cytothreshold=now-t0;
t0=now; 

%% segment cell using nuc and finalize three labels matriceis
clear MAPK
cyto_ringstrel=arg.cyto_ringstrel; 
parfor i=1:numel(T)
    bw = BW(:,:,i) & imdilate(NucLabels(:,:,i)>0,cyto_ringstrel);
    CellLabels(:,:,i) = segmentUsingSeeds(bw,NucLabels(:,:,i)); 
    cyt = CellLabels(:,:,i); 
    cyt(NucLabels(:,:,i)>0)=0; 
    CytoLabels(:,:,i)=cyt; 
    NucLabels(:,:,i) = imerode(NucLabels(:,:,i),strel('disk',nuc_smooth1)); 
end

Timing.cytosegmentation=now-t0;
t0=now; 

%% create the Lbl and populate it
Lbl = CellLabel;
Lbl.saveToFile=true; 
Lbl.pth=MD.pth; 
Lbl.posname = well; 
Lbl.Reg = Reg; 
for i=1:numel(T)
    %% add to Lbl
    addLbl(Lbl,CellLabels(:,:,i),'base',T(i),'relabel','none');
    addLbl(Lbl,CytoLabels(:,:,i),'cyto',T(i),'relabel','none');
    addLbl(Lbl,NucLabels(:,:,i),'nuc',T(i),'relabel','none');
end
   
Timing.createlabelobject=now-t0;

