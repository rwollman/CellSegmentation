function Lbl = segmenetEKARev(MD,well,varargin)
%  Lbl = segmenetEKARev(MD,well,varargin) will perform a segmentation of
%  cells expressing YFP in the cytoplasm and deepblue in the nucleus. Gets
%  as inputs the Metadata object MD, the well to work on and a whole bunch
%  of parameters to define the image analysis routine with comments below.

%% input parameters
arg.verbose = true;
t0=now;

% timefunc provides a mechanism to chose only a subset of images to work on
% based on time timefunc should be a function of time
arg.timefunc = @(t) true(size(t));

arg.projectnucandcyto = false;

arg.positiontype = 'Position'; 

% parameters for nuclei detection
arg.nuc_channel = 'DeepBlue'; 
arg.nuc_erode = strel('disk',3); % initial erosion to enhance nuclei centers
arg.nuc_smooth = fspecial('gauss',7,5); % filtering to smooth it out
arg.nuc_suppress = 0.01; % supression of small peaks - units are in a [0 1] space
arg.nuc_minarea = 30; % smaller then this its not a nuclei

% parameters for cell detection
arg.cyto_dilate = strel('disk',35);
arg.cyto_minarea = 100; % will remove cells with cytoplasm area smaller then this
arg.cyto_maxarea = 2000; % will remove cells with cytoplasm area bigger then this
arg.cyto_minyfp = 0.01; % mimimal yfp intensity in the cytoplasm (absolute units, euqal to 2^16 of ~650)
arg.cyto_distfromedge=100; % remove cells that are in the edge of the image
arg.cyto_thresh = @median; % how to determine the intracellular threshold

arg.replicatenuc = false; 

arg = parseVarargin(varargin,arg);

%% Create the CellLabel object
Lbl = CellLabel;
% 
% %% get timepoints for the Label matrices
% T = MD.getSpecificMetadata('TimestampFrame','Channel','Yellow','Position',well,'timefunc',arg.timefunc);
% T = cat(1,T{:});

%% read data
% read Hoecht images and find the corresponding yellow image using the
% frame Timestamp;
T = MD.getSpecificMetadata('TimestampFrame','Channel',arg.nuc_channel,arg.positiontype,well,'timefunc',arg.timefunc);
T=cat(1,T{:});
nuc = stkread(MD,'Position',well,'Channel',arg.nuc_channel,'TimestampFrame',T);
yfp = stkread(MD,'Position',well,'Channel','Yellow','TimestampFrame',T);

if arg.replicatenuc && size(nuc,3)==1
    yfp = stkread(MD,arg.positiontype,well,'Channel','Yellow');
    nuc = repmat(nuc,[1 1 size(yfp,3)]); %changed rempat to repmat
    T = MD.getSpecificMetadata('TimestampFrame','Channel','Yellow',arg.positiontype,well,'timefunc',arg.timefunc);
    T=cat(1,T{:});
end

[~,ordr]=sort(T); 
nuc=nuc(:,:,ordr); 
yfp=yfp(:,:,ordr); 

if arg.projectnucandcyto
    nuc = mean(nuc,3);
    yfp = mean(yfp,3);
end

CellLabels = cell(size(yfp,3),1);
CytoLabels = cell(size(yfp,3),1);
NucLabels = cell(size(yfp,3),1);

arg.verbose && fprintf('Finishd reading YFP & Hoescht T=%s\n',datestr(now-t0,13));  %#ok<*VUNUS>

assert(size(yfp,3)==size(nuc,3),'Sizes of YFP and DeepBlue should be the same!');

parfor i=1:size(yfp,3)
    %% segment nucleus
    % overall strategy -
    % perform a few operations to transform the nuclei into small hills
    % then find the local maxima to create a seed per each nucleur. Then
    % use watershed to segment the nuclei based on these seeds.
    
    nucprj = nuc(:,:,i);
    msk = nucprj>prctile(nucprj(:),5);
    nucbw = optThreshold(nucprj,'method','otsu','msk',logical(msk),'transform','log');
    
    % identify peaks in the nuc
    
    % subtract bacgkround using a mask to avoid corner issues
    nucprj = backgroundSubtraction(nucprj,'msk',logical(msk));
    
    % erode to emphasize the peaks
    nucpeaks = imerode(nucprj,arg.nuc_erode);
    
    % smooths will get read of all the high frequency information and will
    % leabe a single peak per nuclei. We smooth with an input parametner
    % arg.nuc_smooth
    nucpeaks = imfilter(nucpeaks,arg.nuc_smooth);
    
    % next step is supression by intensity, to move to relative units, I am
    % rescaleing the image such that [0 1] are the 1% and 99% grayscale
    % values
    nucpeaks = mat2gray(nucpeaks,double(prctile(nucprj(nucbw),[1 99])));
    
    % imhmax supresses peaks below a specific height
    nucpeaks = imhmax(nucpeaks,arg.nuc_suppress);
    
    % at this point in the code each nuclei should have a unique sinlgle
    % peak, so we will use a regional max operation to find it. The &nucbw
    % is just to make sure that the identified peaks are within area of
    % hoecht and not some random peak in the background.
    nucpeaks = imregionalmax(nucpeaks) & nucbw;
    
    % next function is a custom segmentation function Roy wrote. It uses
    % combination of bwdistgeodesic and watershed to segment and then
    % removes objects that are too small and deals with the numbering
    nuclbl = segmentUsingSeeds(nucbw,bwlabel(nucpeaks),...
        'method','watershed','mincellarea',arg.nuc_minarea);
    
    arg.verbose && fprintf('%g/%g - Finishd Segmenting nuclei T=%s\n',i,size(yfp,3),datestr(now-t0,13));  %#ok<*VUNUS>
    
    %% segment cytoplasm
    
    % overall strategy - 1. segment / 2. watershed.
    % however, we can't just segment based on forground/background with
    % optThresh since the entire image is foreground. So we relay on cell
    % intrinsic thresholding based on combination of proximity to nuclei
    % and yfp intensity that is higher then the median of the nucleus.
    
    % create an area that surrounds cells with imclose, since nuclei are
    % close to each other this will create a continus true in areas with
    % cells and false in areas without cells.
    possiblecellbw = imclose(nuclbl>0,arg.cyto_dilate);
    
    % now use watershed to divide the foreground area into possible cells
    posscytolbl = segmentUsingSeeds(possiblecellbw,nuclbl,'method','watershed');
    
    % next few steps are use the yfp intensity
    yfpprj = yfp(:,:,i);
    
    % get pixel values for all labels
    PxlIdx_posscyto = regionprops(posscytolbl,'PixelIdxList');
    PxlIdx_posscyto = {PxlIdx_posscyto.PixelIdxList};
    PxlIdx_nuc = regionprops(nuclbl,'PixelIdxList');
    PxlIdx_nuc = {PxlIdx_nuc.PixelIdxList};
    
    % init a false matrix for real cyto
    cytobw = false(size(yfpprj));
    % for each possible cells, find out the median yfp intensity in the
    % nucleus (cyto_thresh is a function to determine threshold, defaults
    % to median) then threshold the possible pixels that might be a cell
    % based on that threshold. This effectivly defines a cells as the
    % pixels that "belog" to a nuclei based on proximity that are also
    % higher in intensity then the median of the nucleus.
    for j=1:numel(PxlIdx_posscyto)
        thrsh = arg.cyto_thresh(yfpprj(PxlIdx_nuc{j}));
        ix = yfpprj(PxlIdx_posscyto{j})>thrsh;
        cytobw(PxlIdx_posscyto{j}(ix)) = true;
    end
    
    % now that we have the trye
    cytolbl = segmentUsingSeeds(cytobw,nuclbl,'method','watershed');
    
    % remove nuclei from cytoplasm
    cytolbl(imerode(nuclbl>0,strel('disk',1)))=0;
    PxlIdx_cyto = regionprops(cytolbl,'PixelIdxList');
    PxlIdx_cyto = {PxlIdx_cyto.PixelIdxList};
    
    %% remove cells by critria
    
    % find for each cell its area, centroid and yfp intensity
    prps = regionprops(cytolbl,yfpprj,{'Area','Centroid','MeanIntensity'});
    Area = [prps.Area];
    YFP = [prps.MeanIntensity];
    
    % find out for each cell its deistance from edge
    bnd=bwboundaries(msk);
    bnd=bnd{1}(:,[2 1])';
    cntr = cat(1,prps.Centroid)';
    DistFromEdge = min(distance(cntr,bnd),[],2);
    
    % combine all criteral to one true/false index & remove cells
    ix = Area>arg.cyto_minarea  & Area<arg.cyto_maxarea  & ...
        YFP > arg.cyto_minyfp & DistFromEdge' > arg.cyto_distfromedge;
    
    PxlIdx_cyto(~ix)=[];
    PxlIdx_nuc(~ix)=[];
    
    n1=numel(PxlIdx_nuc);
    n2=numel(PxlIdx_cyto);
    assert(n1==n2,'Different number of nuclei & cytoplasms - check for error');
    
    %% relabel all three regions into the same
    cytolbl = zeros(size(msk));
    nuclbl=zeros(size(msk));
    celllbl = zeros(size(msk));
    for j=1:numel(PxlIdx_cyto)
        cytolbl(PxlIdx_cyto{j})=j;
        nuclbl(PxlIdx_nuc{j})=j;
        celllbl([PxlIdx_cyto{j}; PxlIdx_nuc{j}])=j;
    end
    CellLabels{i}=celllbl;
    CytoLabels{i}=cytolbl;
    NucLabels{i}=nuclbl;
    
    %% verbose
    arg.verbose && fprintf('%g/%g - Finishd Segmenting cells & relabeling T=%s\n',i,size(yfp,3),datestr(now-t0,13));  %#ok<*VUNUS>
    
end

for i=1:numel(CellLabels)
    %% add to Lbl
    addLbl(Lbl,CellLabels{i},'base',T(i),'relabel','nearest');
    addLbl(Lbl,CytoLabels{i},'cyto',T(i),'relabel','nearest');
    addLbl(Lbl,NucLabels{i},'nuc',T(i),'relabel','nearest');
    
end



