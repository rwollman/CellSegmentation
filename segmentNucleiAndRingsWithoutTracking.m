function [Lbl,NucLabels,CellLabels] = segmentNucleiAndRingsWithoutTracking(MD,well,varargin)

arg.nuc_channel = 'DeepBlue'; 
arg.nuc_smooth1 = fspecial('gauss',9,5);
arg.nuc_smooth2 = fspecial('gauss',15,9);
arg.cyto_channels = {'Cyan','eGFP','Red'}; 
arg.positiontype = 'Position'; 
arg.register = []; % optional registration object
arg.timefunc = @(t) true(size(t));
arg.cyto_ringstrel = strel('disk',15); 
arg = parseVarargin(varargin,arg); 
    
%% Create the CellLabel object
Lbl = CellLabel;

%% get timepoints for the Label matrices
T = MD.getSpecificMetadata('TimestampFrame','Channel',arg.nuc_channel,arg.positiontype,well,'timefunc',arg.timefunc);
T = cat(1,T{:});
nuc = stkread(MD,'TimestampFrame',T,'Channel',arg.nuc_channel); 
[T,ordr]= sort(T);
nuc=nuc(:,:,ordr); 

%% register if requested or if Registration object was passed as an input arguemtn
if islogical(arg.register) && arg.register
    [nuc,Tforms] = registerStack(nuc);
    Reg = Registration(T,Tforms);
    arg.register = Reg; 
elseif ~isempty(arg.register) && isa(arg.register,'Registration')
    Reg = arg.register; 
    nuc = Reg.register(nuc,T); 
else
    Reg = []; 
end

%% first subtrack background for entire stack
% subtract bacgkround using a mask to avoid corner issues
nucprj = mean(nuc,3);
msk = nucprj>prctile(nucprj(:),5);
msk = imfill(msk,'holes');
msk  = bwareaopen(msk,10000);

nuc = backgroundSubtraction(nuc,'msk',logical(msk),'smoothstk',false);
NucLabels = cell(size(nuc,3),1);
CellLabels = cell(size(nuc,3),1);
 
%% for each well, 
parfor i=1:numel(T)
    %% segment all nuclei
    nucfilt = imfilter(nuc(:,:,i),arg.nuc_smooth1) ;
    bw = optThreshold(nucfilt,'msk',msk,'method','localotsu','transform','none');
    pks = imregionalmax(imfilter(nucfilt,arg.nuc_smooth2));
    bw  = bw &~bwareaopen(bw,1000);
    NucLabels{i} = segmentUsingSeeds(bw,bwlabel(pks & bw));
    %% segment cytoplasm
    nuclbl = NucLabels{i};
    cyto = nan(nnz(msk),3); 
    for j=1:numel(arg.cyto_channels)
        c = imread(MD,'TimestampFrame',T(i),'Channel',arg.cyto_channels(j));
        c = backgroundSubtraction(c,'msk',msk); 
        c(~msk)=[]; 
        cyto(:,j)=c; 
    end
    gm = gmdistribution.fit(cyto,2);
    bwfp = zeros(size(nuclbl));
    bwfp(msk)=gm.cluster(cyto)-1; 
    % need to verify which component and finish the clustering... 
    bw = bwfp & imdilate(nuclbl>0,arg.cyto_ringstrel);
    CellLabels{i} = segmentUsingSeeds(bw,nuclbl); 
end
    
