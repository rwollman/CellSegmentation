function Lbl = segmentNucleiOnly(MD,well,varargin)

%% define analysis parameters
% parameters for nuclei detection
arg.nuc_erode = strel('disk',3); % initial erosion to enhance nuclei centers
arg.nuc_smooth = fspecial('gauss',7,5); % filtering to smooth it out
arg.nuc_suppress = 0.01; % supression of small peaks - units are in a [0 1] space
arg.nuc_minarea = 30; % smaller then this its not a nuclei
arg.mindistancefromedge = 150;

arg.positiontype = 'Position'; 

arg = parseVarargin(varargin,arg); 
    
%% read the Hoescht stack
nuc = stkread(MD,arg.positiontype,well,'Channel','DeepBlue');

%% Create the CellLabel object
Lbl = CellLabel;

%% get timepoints for the Label matrices
T = MD.getSpecificMetadata('TimestampFrame','Channel','DeepBlue',arg.positiontype,well);
T = cat(1,T{:});
T = sort(T);

%% first subtrack background for entire stack
% subtract bacgkround using a mask to avoid corner issues
nucprj = mean(nuc,3);
msk = nucprj>prctile(nucprj(:),5);
msk = imfill(msk,'holes');
msk  = bwareaopen(msk,10000);

bnd=bwboundaries(msk);
bnd=bnd{1}(:,[2 1]);

nuc = backgroundSubtraction(nuc,'msk',logical(msk));
NucLabels = cell(size(nuc,3),1);

%% for each well, segment all nucleri
parfor i=1:size(nuc,3)
    %%
    %% segment nucleus
    % overall strategy -
    % perform a few operations to transform the nuclei into small hills
    % then find the local maxima to create a seed per each nucleur. Then
    % use watershed to segment the nuclei based on these seeds.
    nucbw = optThreshold(nuc(:,:,i),'method','otsu','msk',logical(msk),'transform','log');

    % identify peaks in the nuc

    % erode to emphasize the peaks
    nucpeaks = imerode(nuc(:,:,i),arg.nuc_erode);

    % smooths will get read of all the high frequency information and will
    % leabe a single peak per nuclei. We smooth with an input parametner
    % arg.nuc_smooth
    nucpeaks = imfilter(nucpeaks,arg.nuc_smooth);
    
    % next step is supression by intensity, to move to relative units, I am
    % rescaleing the image such that [0 1] are the 1% and 99% grayscale
    % values
    nucslice = nuc(:,:,i);
    nucpeaks = mat2gray(nucpeaks,double(prctile(nucslice(nucbw),[1 99])));
    
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
    
    % remove cells too close to boundaries
    prps = regionprops(nuclbl,'Centroid');
    cntr = cat(1,prps.Centroid)';
    DistFromEdge = min(distance(cntr,bnd'),[],2);
    nucbw = ismember(nuclbl,find(DistFromEdge>arg.mindistancefromedge));
    nucbw = bwareaopen(nucbw,arg.nuc_minarea);
    nuclbl = bwlabel(nucbw);

    NucLabels{i} = nuclbl; 
    
end

%% add label to CellLabel object
for i=1:numel(NucLabels)
    addLbl(Lbl,NucLabels{i} ,'base',T(i),'relabel','nearest');
end