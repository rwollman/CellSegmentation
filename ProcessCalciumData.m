%{
 This function processes the data read the the directory and creates a
 results object 
 
Input: pth- the directory of the raw data 
       pos(optional): a cel array taht contains the specific positions to be analyzed. If this
       argument is empty then all the wells will be analyzed 

Output:
       R - the results object

%}
function [ R ] = ProcessCalciumData( pth,varargin )

    %% define analysis parameters
    % parameters for nuclei detection
    arg.nuc_erode = strel('disk',3); % initial erosion to enhance nuclei centers
    arg.nuc_smooth = fspecial('gauss',7,5); % filtering to smooth it out
    arg.nuc_suppress = 0.01; % supression of small peaks - units are in a [0 1] space
    arg.nuc_minarea = 30; % smaller then this its not a nuclei
    arg.mindistancefromedge = 150;

    %%
    md = Metadata(pth);
    R = MultiPositionSingleCellResults(pth);
    if isempty(varargin)
        R.PosNames = unique(md,'Position');
    else
        R.PosNames = varargin{:};
    end
   
    %% main analysis loop
    for j=1:numel(R.PosNames)

        %% read the Hoescht stack
        nuc = stkread(md,'Position',R.PosNames{j},'Channel','DeepBlue');
        %% Create the CellLabel object
        Lbl = CellLabel;

        %% get timepoints for the Label matrices
        T = md.getSpecificMetadata('TimestampFrame','Channel','DeepBlue','Position',R.PosNames{j});
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

        nuc = backgroundSubtraction(nuc,'msk',logical(msk),'timefilter',sum(fspecial('gauss',50,25)));

        %% register all stack
        [nuc,tforms] = registerStack(nuc,'crop',[680 680 680 680]);
        %[nuc,tforms] = registerStack(nuc,'crop',[1000 1000 1000 1000]);
        % add the Hoescht stack to the results object
        add(R,'nucStack',nuc, R.PosNames{j});
        %% for each well, segment all nucleri
        for i=1:size(nuc,3)
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



            %% add label to CellLabel object
            addLbl(Lbl,nuclbl,'base',T(i),'relabel','nearest');
        end


        %% read Yellow
        ca = stkread(md,'Position',R.PosNames{j},'Channel','Yellow');
        ca1 = ca(:,:,1);
        msk = ca1>prctile(ca1(:),5);

        %% background subtract
        ca = backgroundSubtraction(ca,'msk',logical(msk),'timefilter',sum(fspecial('gauss',50,25)));

        %% get the Time vector
        Tca = md.getSpecificMetadata('TimestampFrame','Channel','Yellow','Position',R.PosNames{j});
        Tca = cat(1,Tca{:});
        [Tca ,ordr]= sort(Tca);

        %% register Ca stack
        ref2d = imref2d([size(ca,1) size(ca,2)]);
        for i=1:size(Tca)
            % find the closest transformation (as calculated by nuclei) and use
            % it to adjust the ca image
            [~,ix] = min(abs(T-Tca(i)));

            % spatially transform (warp)
            ca(:,:,i) = imwarp(ca(:,:,i),tforms(ix),'OutputView',ref2d);
        end
        %add the Ca stack to the results object
        %add(R,'CaStack',ca,R.PosNames{j});
        %% measure mean intensity
        ca = ca(:,:,ordr);
        Ca = meanIntensityPerLabel(Lbl,ca,Tca,'func','mean','type','base');

        %% add to Results object
        setLbl(R,Lbl,R.PosNames{j});
        addTimeSeries(R,'Ca',Ca,Tca,R.PosNames{j});

        %% add Ca2+ stack
        % resize ca stack into 0.25 of its size with stkfun
        % add using R.addData(ca,R.PosNames{j})

        %{


        %make movies
        switch R.PosNames{j}

            case 'B07'
                moviename = 'ATP_Ca_Response_10um_1st';
            case 'B03'
                moviename = 'ATP_Ca_Response_10um_2nd';
            case 'C03'
                moviename = 'ATP_Ca_Response_10um_3rd';
            case 'B08'
                moviename = 'ATP_Ca_Response_1um_1st';
            case 'B04'
                moviename = 'ATP_Ca_Response_1um_2nd';
            case 'C04'
                moviename = 'ATP_Ca_Response_1um_3rd';
            case 'B09'
                moviename = 'ATP_Ca_Response_0.1um_1st';
            case 'B05'
                moviename = 'ATP_Ca_Response_0.1um_2nd';
            case 'C05'
                moviename = 'ATP_Ca_Response_0.1um_3rd';
            case 'B10'
                moviename = 'ATP_Ca_Response_Control';
        end

        R.stack2movie(ca,moviename,'colormap',gray(256));
        %}
    end

    % %%
    %  ix = unidrnd(size(Ca,2));
    %  plot(Tca,Ca(:,ix),'x-'),
    %  ylim([0.004 0.01]),
    %  title(sprintf('X: %g Y: %g',XY(ix,1),XY(ix,2)))


 
end

