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
arg.shrinkmsk = strel('disk',50);
arg.ring_spacer = strel('disk',1); 
arg.ring_width = strel('disk',4);
arg.sz =  [2048  2064]; 
arg.track = 'none'; 
arg.verbose = false; 

arg = parseVarargin(varargin,arg); 
    
t0=now; 
t00=now; 

%% get timepoints for the Label matrices
T = MD.getSpecificMetadata('TimestampFrame','Channel',arg.nuc_channel,arg.positiontype,well,'timefunc',arg.timefunc);
T = cat(1,T{:});
nuc = stkread(MD,'TimestampFrame',T,'Channel',arg.nuc_channel); 
[T,ordr]= sort(T);
nuc=nuc(:,:,ordr); 

Timing.readnucimages=now-t0;
t0=now; 

if arg.register==1
    arg.register=true; 
end
if arg.register==0
    arg.register=false; 
end

%% register if requested or if Registration object was passed as an input arguemtn
if islogical(arg.register) && arg.register
    
    %% find the frames where potential shifts occured (between acq)
    tbl=MD.tabulate('acq','Position',well,'Channel',arg.nuc_channel,'timefunc',arg.timefunc); 
    tbl=cat(1,tbl{:,2}); 
    possibleShiftingFrames=cumsum(tbl(1:end-1))+1; 
    [nuc,Tforms] = registerStack(nuc,'crop',[400 400 1200 1200],'method','xcorr','reference',1,'maxdisp',100,'onlyspecificframes',possibleShiftingFrames);
    Reg = Registration(T,Tforms);
    arg.register = Reg; 
elseif ~isempty(arg.register) && isa(arg.register,'Registration')
    Reg = arg.register; 
    nuc = Reg.register(nuc,T); 
else
    Reg = []; 
end

Timing.registration=now-t0;
arg.verbose && fprintf('Finished registration %s\n',datestr(now-t00,13)); 

%% first subtrack background for entire stack
% subtract bacgkround using a mask to avoid corner issues
nucprj = mean(nuc,3);
msk = nucprj>prctile(nucprj(:),5);
msk = imfill(msk,'holes');
msk  = bwareaopen(msk,10000);
if ~isempty(arg.shrinkmsk)
    msk = imerode(msk,arg.shrinkmsk); 
end
nuc = backgroundSubtraction(nuc,'msk',logical(msk),'smoothstk',false);


NucLabels = zeros([arg.sz size(nuc,3)],'uint16');
CellLabels = zeros([arg.sz size(nuc,3)],'uint16');
CytoLabels = zeros([arg.sz size(nuc,3)],'uint16');
RingLabels = zeros([arg.sz size(nuc,3)],'uint16');

Timing.backgroundsubtraction=now-t0;
t0=now; 
arg.verbose && fprintf('Finished background subtraction %s\n',datestr(now-t00,13)); 

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
arg.verbose && fprintf('Finished nuclei segmentation %s\n',datestr(now-t00,13)); 

%% create a whole cell binary mask
BW = true(size(nuc)); 
for j=1:numel(arg.cyto_channels)
    MAPK = stkread(MD,'TimestampFrame',T,'Channel',arg.cyto_channels{j});
    if ~isempty(Reg)
        MAPK = Reg.register(MAPK,T); 
    end
    MAPK = imfilter(MAPK,fspecial('gauss',5,3)); 
    cyto_thresholdmethod=arg.cyto_thresholdmethod; 
    cyto_transform=arg.cyto_transform; 
    BWtmp=false(size(nuc)); 
    parfor i=1:numel(T)
        %% segment cytoplasm
        BW(:,:,i) = BW(:,:,i) & optThreshold(MAPK(:,:,i),'msk',msk,'method',arg.cyto_thresholdmethod,'transform',arg.cyto_transform);
        BWtmp(:,:,i) = optThreshold(MAPK(:,:,i),'msk',msk,'method',cyto_thresholdmethod,'transform',cyto_transform);
    end
    BW=BW | BWtmp; 
end

Timing.cytothreshold=now-t0;
t0=now; 
arg.verbose && fprintf('Finished cyto thresholding %s\n',datestr(now-t00,13)); 

%% segment cell using nuc and finalize three labels matriceis
clear MAPK
cyto_ringstrel=arg.cyto_ringstrel; 
ring_spacer = arg.ring_spacer; 
ring_width=arg.ring_width; 

parfor i=1:numel(T)
    
    %% limit cell signal to some positive distance from cyto
    nucmskdil =  imdilate(NucLabels(:,:,i)>0,cyto_ringstrel); 
    bw = BW(:,:,i) & nucmskdil;
        
    %% nuc must be where there is a cell signal, i.e. bw is true
    nuclbl = NucLabels(:,:,i); 
    nuclbl(~bw)=0; 
    nuclbl = bwlabel(nuclbl>0); 
    
    celllbl = segmentUsingSeeds(bw,nuclbl); 
    cytlbl = celllbl; 
    cytlbl(nuclbl>0)=0; 
    
      
    %% some bookeeping to make sure the numbers of all three match
    % get rid of nuclei without cytoplasm
    nuclbl(ismember(nuclbl,setdiff(nuclbl,cytlbl)))=0; 
    nuclbl = bwlabel(nuclbl>0); 
    
    % redo segmentation but now with only "correct" nuclei
    celllbl = segmentUsingSeeds(bw,nuclbl); 
    cytlbl = celllbl; 
    cytlbl(nuclbl>0)=0; 
    
    % create a ring
    nuclbl2 = imdilate(nuclbl,ring_spacer);
    rnglbl=imdilate(nuclbl2,ring_width);
    rnglbl(nuclbl2>0)=0; % remove nuclei
    rnglbl(cytlbl==0)=0; % remove ring pixesl outside of cytoplasm
        
    
    % assign everything
    CellLabels(:,:,i)=celllbl;
    CytoLabels(:,:,i)=cytlbl; 
    NucLabels(:,:,i) = nuclbl; 
    RingLabels(:,:,i) = rnglbl; 
end

Timing.cytosegmentation=now-t0;
t0=now; 
arg.verbose && fprintf('Finished cyto segmentation %s\n',datestr(now-t00,13)); 

%% create the Lbl and populate it
Lbl = CellLabel;
Lbl.saveToFile=true; 
Lbl.pth=MD.pth; 
Lbl.posname = well; 
Lbl.Reg = Reg; 
addLbl(Lbl,CellLabels,'base',T,'relabel',arg.track);
for i=1:numel(T)
    addLbl(Lbl,CytoLabels(:,:,i),'cyto',T(i),'relabel',arg.track);
    addLbl(Lbl,NucLabels(:,:,i),'nuc',T(i),'relabel',arg.track);
    addLbl(Lbl,RingLabels(:,:,i),'ring',T(i),'relabel',arg.track);
end
   
arg.verbose && fprintf('Finished creating cell label %s\n',datestr(now-t00,13)); 
Timing.createlabelobject=now-t0;

