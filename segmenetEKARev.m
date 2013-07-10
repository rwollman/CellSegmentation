function Lbl = segmenetEKARev(MD,well,varargin)

%% input parameters
arg.verbose = true; 
t0=now; 

arg.measurefret = true; 
arg.measurexrhod = true; 

arg.timefunc = @(t) true(size(t)); 

% parameters for nuclei detection
arg.nuc_smooth = fspecial('gauss',7,5); 
arg.nuc_suppress = 0.01; 
arg.nuc_erode = strel('disk',3); 
arg.nuc_minarea = 30; 

% parameters for cell detection
arg.cyto_dilate = strel('disk',35); 
arg.cyto_minarea = 100; 
arg.cyto_maxarea = 2000; 
arg.cyto_minyfp = 0.01; 
arg.cyto_distfromedge=100; 
arg.cyto_thresh = @median; 

arg.redo = false; 

arg = parseVarargin(varargin,arg); 

%% Create the CellLabel object
Lbl = CellLabel; 

%% get timepoints for the Label matrices
T = MD.getSpecificMetadata('TimestampFrame','Channel','Yellow','Position',well,'timefunc',arg.timefunc);
T = cat(1,T{:});

%% read data
yfp = stkread(MD,'Position',well,'Channel','Yellow','timefunc',arg.timefunc);
nuc = stkread(MD,'Position',well,'Channel','DeepBlue','timefunc',arg.timefunc);

arg.verbose && fprintf('Finishd reading YFP & Hoescht T=%s\n',datestr(now-t0,13));  %#ok<*VUNUS>

assert(size(yfp,3)==size(nuc,3),'Sizes of YFP and DeepBlue should be the same!'); 

for i=1:size(yfp,3) 

    %% segment nucleus
    nucprj = nuc(:,:,i);
    msk = nucprj>prctile(nucprj(:),5);
    nucbw = optThreshold(nucprj,'method','otsu','msk',logical(msk),'transform','log');
    
    % identify peaks in the nuc
    nucprj = backgroundSubtraction(nucprj,'msk',logical(msk));
    nucpeaks = imerode(nucprj,arg.nuc_erode);
    nucpeaks = imfilter(nucpeaks,arg.nuc_smooth);
    nucpeaks = mat2gray(nucpeaks,double(prctile(nucprj(nucbw),[1 99])));
    nucpeaks = imhmax(nucpeaks,arg.nuc_suppress);
    nucpeaks = imregionalmax(nucpeaks) & nucbw;
    nuclbl = segmentUsingSeeds(nucbw,bwlabel(nucpeaks),...
        'method','watershed','mincellarea',arg.nuc_minarea);
    
    arg.verbose && fprintf('%g/%g - Finishd Segmenting nuclei T=%s\n',i,size(yfp,3),datestr(now-t0,13));  %#ok<*VUNUS>
                       
    %% segment cytoplasm
    possiblecellbw = imclose(nuclbl>0,arg.cyto_dilate);
    posscytolbl = segmentUsingSeeds(possiblecellbw,nuclbl,'method','watershed');
    yfpprj = yfp(:,:,i);
    
    PxlIdx_posscyto = regionprops(posscytolbl,'PixelIdxList');
    PxlIdx_posscyto = {PxlIdx_posscyto.PixelIdxList};
    PxlIdx_nuc = regionprops(nuclbl,'PixelIdxList');
    PxlIdx_nuc = {PxlIdx_nuc.PixelIdxList};
    cytobw = false(size(yfpprj));
    for j=1:numel(PxlIdx_posscyto)
        thrsh = arg.cyto_thresh(yfpprj(PxlIdx_nuc{j}));
        ix = yfpprj(PxlIdx_posscyto{j})>thrsh;
        cytobw(PxlIdx_posscyto{j}(ix)) = true;
    end
    cytolbl = segmentUsingSeeds(cytobw,nuclbl,'method','watershed');
    cytolbl(imerode(nuclbl>0,strel('disk',1)))=0;
    PxlIdx_cyto = regionprops(cytolbl,'PixelIdxList');
    PxlIdx_cyto = {PxlIdx_cyto.PixelIdxList};

    %% remove cells by critria

    prps = regionprops(cytolbl,yfpprj,{'Area','Centroid','MeanIntensity'});
    Area = [prps.Area];
    YFP = [prps.MeanIntensity];
    
    bnd=bwboundaries(msk);
    bnd=bnd{1}(:,[2 1])';
    cntr = cat(1,prps.Centroid)';
    DistFromEdge = min(distance(cntr,bnd),[],2);
    
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

    %% add to Lbl
    addLbl(Lbl,celllbl,'base',T(i),'relabel','nearest');
    addLbl(Lbl,cytolbl,'cyto',T(i),'relabel','nearest'); 
    addLbl(Lbl,nuclbl,'nuc',T(i),'relabel','nearest');
    
    %% verbose
    arg.verbose && fprintf('%g/%g - Finishd Segmenting cells & relabeling T=%s\n',i,size(yfp,3),datestr(now-t0,13));  %#ok<*VUNUS> 
end


   

