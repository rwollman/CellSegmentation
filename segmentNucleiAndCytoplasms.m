function Lbl = segmentNucleiAndCytoplasms(MD,well,varargin)

%% input parameters
arg.verbose = true;
t0=now;

% timefunc provides a mechanism to chose only a subset of images to work on
% based on time timefunc should be a function of time
arg.timefunc = @(t) true(size(t));

arg.projectnucandcyto = false;

% parameters for nuclei detection - mostly critical for under/over
% segmentation
arg.nuc_erode = strel('disk',3); % initial erosion to enhance nuclei centers
arg.nuc_smooth = fspecial('gauss',7,5); % filtering to smooth it out
arg.nuc_suppress = 0.01; % supression of small peaks - units are in a [0 1] space
arg.nuc_minarea = 30; % smaller then this its not a nuclei

arg.nuc_stretch = [1 99]; 
arg.nuc_channel = 'DeepBlue'; 

arg.project = false; 

arg.mindistancefromedge =0;

% parameters for cell detection
arg.cyto_channel = 'Yellow'; 

% chose algorithm to use to threshold the cyto
arg.cyto_threshalgo='nucleibased'; % other option is {'optthresh'};

% relevant for the nucleibased algorithm 
arg.cyto_dilate = strel('disk',35);
arg.cyto_thresh = @median; % how to determine the intracellular threshold

% relevant for the optthresh algorithm
arg.cyto_threshmethod='otsu'; 
arg.cyto_threshtransform='log'; 

arg.cyto_minarea = 100; % will remove cells with cytoplasm area smaller then this
arg.cyto_maxarea = 2000; % will remove cells with cytoplasm area bigger then this
arg.cyto_minyfp = 0; % mimimal yfp intensity in the cytoplasm (absolute units, euqal to 2^16 of ~650)
arg.cyto_distfromedge=100; % remove cells that are in the edge of the image

arg.positiontype = 'Position'; 

arg.reg_channel = 'CyanToYellow'; 
arg.register = []; 

arg.shrinkmsk = strel('disk',50);

arg.track_method='none'; %nearest %none

arg = parseVarargin(varargin,arg);

%% Create the CellLabel object

%% Create and populate the registration object
if ~arg.register
    Reg=false; 
elseif isa(arg.register,'Registration')
    Reg = arg.register; 
else
    Reg = Registration;
    Treg = MD.getSpecificMetadata('TimestampImage',arg.positiontype ,well,'Channel',arg.reg_channel,'timefunc',arg.timefunc);
    Treg = sort(cat(1,Treg{:}));
    Reg.T=Treg;
    stk = stkread(MD,arg.positiontype ,well,'Channel',arg.reg_channel,'timefunc',arg.timefunc);
    [~,Tforms] = registerStack(stk);
    Reg.Tforms = Tforms;
    % add the registration object to the arg struct
    arg.register = Reg;
end



%% call the nuclei segmentation function with updated arg struct
[~,NucLabels,T,msk] = segmentNucleiOnly(MD,well,arg); 
% n = numel(MD.getSpecificMetadata('TimestampImage','Channel',arg.cyto_channel,arg.positiontype,well));



%% read cytoplasm channel and register it with Reg
[yfp,indx] = stkread(MD,'Channel',arg.cyto_channel,arg.positiontype,well);


% get Tyfp for the cytoplasm based on indexes where Tyfp is based on Frame
Tyfp = MD.getSpecificMetadataByIndex('TimestampFrame',indx); 
Tyfp = cat(1,Tyfp{:}); 

if isa(Reg,'Registration') % if Reg returns not false its a 
    yfp = Reg.register(yfp,Tyfp); 
end

%% remove extra timepoints from the nuclabel if there are and duplicate the nuclbl as needed. 
% where we drop any nuclbl where the Tyfp 
% neighbor gets dropped 
NucLabelToKeep = ismember(T,Tyfp); 
NucLabels = NucLabels(NucLabelToKeep); 
T = T(NucLabelToKeep); 

tfYfpWithNuc = ismember(Tyfp,T);
NucLabelsAll = cell(numel(Tyfp),1); 
NucLabelsAll(tfYfpWithNuc)=NucLabels; 
ixMissingNucs = find(cellfun(@isempty,NucLabelsAll)); 
 
for i=1:numel(ixMissingNucs)
    dT = abs(T-Tyfp(ixMissingNucs(i)));  
    [~,mi]=min(dT); 
    NucLabelsAll{ixMissingNucs(i)} = NucLabels{mi}; 
end
NucLabels = NucLabelsAll; 

n = numel(NucLabels); 
CellLabels = cell(n,1);
CytoLabels = cell(n,1);

if numel(T)<numel(Tyfp) 
    T=Tyfp;
end

%% segment cytoplasm
for i=1:numel(NucLabels)
    nuclbl = NucLabels{i}; 
    % overall strategy - 1. segment / 2. watershed.
    % however, we can't just segment based on forground/background with
    % optThresh since the entire image is foreground. So we relay on cell
    % intrinsic thresholding based on combination of proximity to nuclei
    % and yfp intensity that is higher then the median of the nucleus.
     
    % next few steps are use the yfp intensity
    yfpprj = yfp(:,:,i);
    % init a false matrix for real cyto
    switch arg.cyto_threshalgo  
       
        case 'nucleibased'
            % create an area that surrounds cells with imclose, since nuclei are
            % close to each other this will create a continus true in areas with
            % cells and false in areas without cells.
            possiblecellbw = imclose(nuclbl>0,arg.cyto_dilate);
            
            % now use watershed to divide the foreground area into possible cells
            posscytolbl = segmentUsingSeeds(possiblecellbw,nuclbl,'method','watershed');
            
            
            
            % get pixel values for all labels
            PxlIdx_posscyto = regionprops(posscytolbl,'PixelIdxList');
            PxlIdx_posscyto = {PxlIdx_posscyto.PixelIdxList};
            PxlIdx_nuc = regionprops(nuclbl,'PixelIdxList');
            PxlIdx_nuc = {PxlIdx_nuc.PixelIdxList};
            
            
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
        case 'optthresh'
            cytobw=optThreshold(yfpprj,'msk',msk,'method',arg.cyto_threshmethod,'transform',arg.cyto_threshtransform);
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

    %% verbose
    arg.verbose && fprintf('%g/%g - Finishd Segmenting cells & relabeling T=%s\n',i,size(yfp,3),datestr(now-t0,13));  %#ok<*VUNUS>
    
end

%% create the Lbl and populate it
Lbl = CellLabel;
Lbl.Reg = Reg; 
for i=1:numel(CellLabels)
    %% add to Lbl
    addLbl(Lbl,CellLabels{i},'base',T(i),'relabel',arg.track_method);
    addLbl(Lbl,CytoLabels{i},'cyto',T(i),'relabel',arg.track_method);
    addLbl(Lbl,NucLabels{i},'nuc',T(i),'relabel',arg.track_method);
end