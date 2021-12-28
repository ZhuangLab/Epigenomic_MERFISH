%% Measure the MERFISH and cell identities
%% Tian LU
%% 3/1/2017
%% Setup path and parameters
% Define data path
dataPath='';
analysisSavePath = SetFigureSavePath([dataPath,'\FOV\'], ...
     'makeDir', true);
RawPath = ''; 
DAPIPath=''; 
% Useful data structure for spotfinding
all=readtable('');

[~,~,codebook]=xlsread('');
code=cell2mat(codebook(1:366,5:28));
Name=codebook(1:366,1);
Gene=codebook(1:366,4);
color=all.color;
desiredround=all.imagingRound;
bitnumber=24;
roundnumber=8;
znumber=30;

mapPath=SetFigureSavePath([analysisSavePath,'\loci_localization\'], ...
     'makeDir', true);
decodePath=SetFigureSavePath([analysisSavePath,'\decode_per_cell\'], ...
     'makeDir', true);
%% ------------------------------------------------------------------------
% Start logging
%%-------------------------------------------------------------------------
if ~isempty(mfilename) % Only log if running as a script
    diaryFile = [analysisSavePath 'matlab_output.log']; % File name
    diary off; % Turn off diary if already in use
    if exist(diaryFile) % Delete existing file
        diary off;
        delete(diaryFile);
    end
    diary(diaryFile); % Set diary file
    diary on;

    % Display information
    PageBreak();
    display(['Running: ' mfilename('fullpath') '.m']);
    display(['Started: ' datestr(now)]);
    
    % Archive script
    copyfile( [mfilename('fullpath'),'.m'],[analysisSavePath,mfilename,'.m']);
    display('------------------------------------------------------------------');
    display(['Copied analysis script to ' analysisSavePath,mfilename,'.m']);
    
    % Start script timer
    scriptTimer = tic;
end
%% Correct the illumination of 647, 750, 561,405

if exist([dataPath 'illumination.mat'],'file')
    load([dataPath 'illumination']);
else
    Path=[RawPath ];
    file=['zscan_(?<fov>[0-9]+)_0'];    
    tempFiles = BuildFileStructure(Path, ...
    'fileExt', 'dax', ...
    'regExp', file, ...
    'fieldNames', {'fov','round'}, ...
    'fieldConv', {@str2num});

    display(['Found ' num2str(length(tempFiles)) ' dax files']);

    % Run analysis of all fov
    tempsum=zeros(2048,2048);
    for f=1:length(tempFiles)
        data=ReadDax(tempFiles(f).filePath,'startFrame', znumber/2, 'endFrame', znumber/2+1); 
        tempsum=tempsum+double(data(:,:,1));
        tempsum=tempsum+double(data(:,:,2));       
    end
    
    a=tempsum/length(tempFiles);
    amax=max(max(a));
    ratio750=a/amax;
    

    % Run analysis of all fov
    tempsum=zeros(2048,2048);
    for f=1:length(tempFiles)
        data=ReadDax(tempFiles(f).filePath,'startFrame', znumber*1.5, 'endFrame', znumber*1.5+1); 
        tempsum=tempsum+double(data(:,:,1));
        tempsum=tempsum+double(data(:,:,2));       
    end
    
    a=tempsum/length(tempFiles);
    amax=max(max(a));
    ratio647=a/amax;
    
    tempsum=zeros(2048,2048);
    for f=1:length(tempFiles)
        data=ReadDax(tempFiles(f).filePath,'startFrame', znumber*2.5, 'endFrame', znumber*2.5+1); 
        tempsum=tempsum+double(data(:,:,1));
        tempsum=tempsum+double(data(:,:,2));       
    end
    
    a=tempsum/length(tempFiles);
    amax=max(max(a));
    ratio561=a/amax;
  %%  
    Path=DAPIPath;
    file=['zscan_(?<fov>[0-9]+)_0'];    
    tempFiles = BuildFileStructure(Path, ...
    'fileExt', 'dax', ...
    'regExp', file, ...
    'fieldNames', {'fov','round'}, ...
    'fieldConv', {@str2num});
    % Run analysis of all fov
    tempsum=zeros(2048,2048);
    for f=1:length(tempFiles)
        data=ReadDax(tempFiles(f).filePath,'startFrame', 5, 'endFrame', 5); 
        tempsum=tempsum+double(data(:,:,1));      
    end
    
    a=tempsum/length(tempFiles);
    amax=max(max(a));
    ratio405=a/amax;
 %   

   save([dataPath 'illumination'],'ratio750','ratio647','ratio561','ratio405');

end  
%% Create parallel pool
if isempty(gcp('nocreate'))
    p = parpool(20); % Set this number to control the number of parallel workers.
                     % Polite usage would suggest these maximum values: morgan = 20, cajal = 10
else
    p = gcp;
end


% Try to segment cells on each FOV
nucleusbwThresh1=0.6;
nucleusbwThresh2=0.4;
cytoplasmbwThresh=0.03;

histFcn = @(x) histogram(x, 'Normalization', 'pdf');
load([dataPath 'illumination']);
%%
parfor fovID =1:100
    %%
    if ~exist(strcat(mapPath,['FOV_' num2str(fovID),'data.mat']),'file') 
        Nucleus=ReadDax([DAPIPath,'zscan_',num2str(fovID-1,'%03d'),'_0.dax'], 'startFrame', 5, 'endFrame', 5,'verbose',false);
        if sum(find(Nucleus>300))>100
        Nucleus=imadjust((int16(double(Nucleus)./ratio405)));
%         
%          figure;
%          imshow(imadjust(int16(Nucleus)));

        nucleusImgsi = im2double(Nucleus);
        nucleusImgsi = nucleusImgsi - min(nucleusImgsi(:));
        nucleusImgsi = nucleusImgsi./max(nucleusImgsi(:));
        nucleusbw = (nucleusImgsi>0.55);
        nucleusbw2 = imfill(nucleusbw, 'holes');
        nucleusbw3 = nucleusbw2-bwareaopen(nucleusbw2,5000);
        nucleusbw4=bwareaopen(nucleusbw3,100);
        oldnucleusbw=nucleusbw4;
       
        
        for SPThresh=0.6:0.02:0.9
          nucleusbw = (nucleusImgsi>SPThresh);            
          nucleusbw2 = imfill(nucleusbw, 'holes');
          nucleusbw3 = nucleusbw2-bwareaopen(nucleusbw2,5000);
          nucleusbw4=bwareaopen(nucleusbw3,100)+oldnucleusbw;
          oldnucleusbw=nucleusbw4;
        end
        nucleusbw= imclose(oldnucleusbw, ones(6,6));
        nucleusbw2 = imfill(nucleusbw, 'holes');
        nucleusbw3 = bwareaopen(nucleusbw2,10);
        nucleusbw5 = imopen(nucleusbw3,strel( 'disk', 10 ) );
       % nucleusbw5=imdilate(nucleusbw5,ones(10,10));
        nucleus_perim = bwperim(nucleusbw5);

%         figure;
%         imshow(nucleusbw5);
%         
        %
        cytoplasmImgsi = nucleusImgsi;
        % change cytoplasmImgsi to binary images
        cytoplasmbw = (cytoplasmImgsi >0.5);
        cytoplasmbw1 = imfill(cytoplasmbw,'holes');
        % get rid of the isolated dots and smooth the edge
        cytoplasmbw2 = bwareaopen(cytoplasmbw1, 20);
        % erode the edges and get the eroded image dilated
        cytoplasmbw3 = imdilate(cytoplasmbw2,ones(30,30));
        cytoplasmbw3 = imfill(cytoplasmbw3, 'holes');
%         cytoplasmbw3 = cytoplasmbw3-bwareaopen(cytoplasmbw3,50000);
        % figure; imshowpair(cytoplasmImgsi, cytoplasmbw3, 'montage');
        cytoplasmbw3_perim = bwperim(cytoplasmbw3);
%         
        figure;
        imshow(cytoplasmbw3);
%        
        % find centroid for nucleus
        [labeledImage, numberOfBlobs] = bwlabel(nucleusbw5);%bwlabel(nucleus_perim ~= 0);
        measurements = regionprops(labeledImage, 'Centroid');
        allCentroids = [measurements.Centroid];
        xCentroids = allCentroids(1:2:end);
        yCentroids = allCentroids(2:2:end);

        mask_em = nucleusbw5;
        mask_em = imclose(mask_em, ones(5,5));
        mask_em = imfill(mask_em, 'holes');
        mask_em = bwareaopen(mask_em, 40);

        cytoplasmImgsi_c = imcomplement(imadjust(cytoplasmImgsi));
        I_mod = imimposemin(cytoplasmImgsi_c, ~cytoplasmbw3 | mask_em);
        L = watershed(I_mod);
        cellNum = max(L(:))+1;
        cellCandidates = cell(cellNum,1);
        for celli = 0:cellNum-1
            [cellCandidates{celli+1}(:,1), cellCandidates{celli+1}(:,2)] = ind2sub(size(L), find(L==celli));
        end

        nucleusxCentroids = floor(xCentroids); % convert from float to pixels
        nucleusyCentroids = floor(yCentroids);
        % Count the nucleus number within FOVi
        nucleusCount = length(nucleusyCentroids);
        cellIdx = [];
        nucleusIdx = [];
        nucleusL = labeledImage;%bwlabel(nucleusbw5);
        for celli = 1:cellNum
            % -------------------------------------------------------------------------
            % Keep only identified cell regions with nucleus staining of FOVi in the
            % middle
            % -------------------------------------------------------------------------
            for j = 1:length(nucleusxCentroids)
                distMat =  sqrt((cellCandidates{celli}(:,1) - nucleusyCentroids(j)).^2 + (cellCandidates{celli}(:,2) - nucleusxCentroids(j)).^2);                
                if sum(distMat < 5)
                    cellIdx(end+1) = celli;
                    nucleusIdx(end+1) = j;
                    break;
                end
            end
        end
        cellInfo = cellCandidates(cellIdx);
        nucleusColCentroid = round(nucleusxCentroids(nucleusIdx));
        nucleusRowCentroid = round(nucleusyCentroids(nucleusIdx));
        % reconstruct color map for L to represent cell with >0, and empty regions
        % as 0
        cellMap = zeros(size(L));
        for celli = 1:length(cellInfo)
            cellMap(sub2ind(size(cellMap),cellInfo{celli}(:,1), cellInfo{celli}(:,2))) = celli;
        end
        %     runningTime = toc;
        %     fprintf('Time for selecting out potential cells in the FOVi: %f\n s.', runningTime); 
        % -------------------------------------------------------------------------
        % Plot and keep records of segregated cells
        % -------------------------------------------------------------------------

        cellCount=0;
        cellBoundaryStruct = struct(...
            'cytoplasmBoundary',{},...
            'nucleusBoundary',{},...
            'cellID',[]);
        for celli = 1:length(cellInfo)
            % give cell ID
            cellCount = cellCount + 1;
            cellBoundaryStruct(end+1).cellID = cellCount;
            % cytoplasm boundary
            boundaryPosi =  bwboundaries(cellMap == celli);
            boundaryPosi = boundaryPosi{1};
            cellBoundaryStruct(end).cytoplasmBoundary = boundaryPosi;
             % nucleus boundary
            boundaryPosi = bwboundaries(nucleusL == nucleusIdx(celli));% nucleusL(nucleusRowCentroid(celli), nucleusColCentroid(celli)));
            boundaryPosi = boundaryPosi{1};
            cellBoundaryStruct(end).nucleusBoundary =boundaryPosi;
        end
%         cellBoundaryStructnu=cellBoundaryStruct;
%         %%segment cells
% 
%   figure;
%   imshow(label2rgb(cellMap));


       ValidCell=ones(1,length(cellBoundaryStruct),'int16');
       
  
       %% load raw images
       warning('off','all')
        % start finding spots in 3D
        disp(['Start mapping the spots to cells: FOV',num2str(fovID)]);       
        ratio=zeros(bitnumber,length(cellBoundaryStruct));
        bitinfo=[];       
        
       %map spot based on 561 beads 
        fixed=ReadDax([RawPath,'zscan_',num2str(fovID-1,'%03d'),'_00.dax'], 'startFrame', znumber*3+1, 'endFrame', znumber*4,'verbose',false);
        data=fixed;
        data=flip(data,3);
        h=histFcn(data);
       ccum=zeros(1,h.NumBins);
       for i=1:h.NumBins
           ccum(i)=sum(h.BinCounts(1:i));
       end          
       x=[];
       thre=0.999;
       idx=find(ccum>thre*ccum(end));  
       bw=data>h.BinEdges(idx(1));     
       bw2=imfill(bw,'holes');
       bw3=bwareaopen(bw2,10); 
       Centroid = regionprops3(bw3, 'Centroid');
       cfixed=Centroid.Centroid;
       
        disp(size(cfixed));
       bw4=bw3;
   if sum(ValidCell)>0
       %%
       allspot=cell(1,length(cellBoundaryStruct));
       fixed=[];
       shiftpre=[0,0,0];
       shiftall={};
        for bit=1:bitnumber
            tic;
            shift=[0,0,0];
            switch color(bit)
                case 750
                    a=ReadDax([RawPath,'\zscan_',num2str(fovID-1,'%03d'),'_', num2str(desiredround(bit),'%02d'),'.dax'], 'startFrame', 1, 'endFrame', znumber,'verbose',false);                  
                    a=double(a)./ratio750/100;
                case 650
                    a=ReadDax([RawPath,'\zscan_',num2str(fovID-1,'%03d'),'_', num2str(desiredround(bit),'%02d'),'.dax'], 'startFrame', znumber+1, 'endFrame', znumber*2,'verbose',false);  
                    a=double(a)./ratio647/100;
                    a=flip(a,3);
                case 561
                    a=ReadDax([RawPath,'\zscan_',num2str(fovID-1,'%03d'),'_', num2str(desiredround(bit),'%02d'),'.dax'], 'startFrame', znumber*2+1, 'endFrame', znumber*3,'verbose',false);  
                    a=double(a)./ratio561/100;
            end

            if bit>1 && desiredround(bit)==desiredround(bit-1) 
                   shift=shiftpre;
%                     disp('fov:');
%                     disp(fovID);
            elseif bit>1
                 moving=ReadDax([RawPath,'\zscan_',num2str(fovID-1,'%03d'),'_', num2str(desiredround(bit),'%02d'),'.dax'], 'startFrame', znumber*3+1,'endFrame',znumber*4,'verbose',false);   

                 data=moving;
                 data=flip(data,3);
                   h=histFcn(data);
                   ccum=zeros(1,h.NumBins);
                   for i=1:h.NumBins
                       ccum(i)=sum(h.BinCounts(1:i));
                   end          
                   x=[];
                   thre=0.999;
                   idx=find(ccum>thre*ccum(end));  
                   bw=data>h.BinEdges(idx(1));     
                   bw2=imfill(bw,'holes');
                   bw3=bwareaopen(bw2,10); 
                   Centroid = regionprops3(bw3, 'Centroid');
                   cmove=Centroid.Centroid;
                    moving=[];
                    %%
                    disp('Start warping');
                    cmovet=cmove;
                    
                    for trial=1:10
                         dist=pdist2(cmovet(:,1:2),cfixed(:,1:2));
                        [M,I]=min(dist);
                        [temp,ids]=sort(M);
                        dev=diff(temp);
                        idd=(dev<0.3*median(dev));
                        temp1=cfixed(ids(idd),:);
                        temp2=cmovet(I(ids(idd)),:);
                        co=temp1-temp2;
                    coordi=[median(co(:,1)),median(co(:,2)),median(co(:,3))];    
                    cmovet=[cmovet(:,1)+coordi(1,1),cmovet(:,2)+coordi(1,2),cmovet(:,3)+coordi(1,3)];
                    end

                        shift(1)=mean(cmovet(:,1)-cmove(:,1));
                        shift(2)=mean(cmovet(:,2)-cmove(:,2));
                        shift(3)=mean(cmovet(:,3)-cmove(:,3));
                        shiftpre=shift;
                 end

            disp(['finished warping for FOV ',num2str(fovID), ' bit' num2str(bit)]);
            shiftall{bit}=shift;
            disp(shift)
            %%
%             figure;
%             imshow(imadjust(int16(max(a,[],3))));
%             hold on;
           tempcount=[];
           temptotal=[];
           
                      h=histFcn(a);
                   ccum=zeros(1,h.NumBins);
                   for i=1:h.NumBins
                       ccum(i)=sum(h.BinCounts(1:i));
                   end  
                   temp=max(a,[],3);
%                    figure;
%                                    imshow(int16(temp),[1.2*min(temp(:)),0.5*max(temp(:))]);
%                                    hold on;
           for l=1:length(cellBoundaryStruct)  
                if ValidCell(l)
                  x1=max(min(cellBoundaryStruct(l).cytoplasmBoundary(:,2))-shift(1)-3,1);
                   x2=min(max(cellBoundaryStruct(l).cytoplasmBoundary(:,2))-shift(1)+3,2048);
                   y1=max(min(cellBoundaryStruct(l).cytoplasmBoundary(:,1))-shift(2)-3,1);
                   y2=min(max(cellBoundaryStruct(l).cytoplasmBoundary(:,1))-shift(2)+3,2048);
                               
                   data=a(y1:y2,x1:x2,:);                              
                   idx=find(ccum>0.999999*ccum(end));  
                                   bw=data>h.BinEdges(idx(1));   
                                   %%bw=data>0.9*max(data(:));     
                                   bw2=imfill(bw,'holes');
                                   bw3=bw2-bwareaopen(bw2,1000); 
                                   bw4=bwareaopen(bw3,3);
                                   oldbw=bw4;
                                   thre=0.999;
                %                    if  bit==5 || bit==10 
                %                        thre=0.992;
                %                    end
                                   for threshold=0.99999:-0.0001:thre%bitthreshold(bit)
                                       idx=find(ccum>threshold*ccum(end));  
                                   bw=data>h.BinEdges(idx(1));
                                       %bw=data>threshold*max(data(:));     
                                       bw2=imfill(bw,'holes');
                                       bw3=bw2-bwareaopen(bw2,1000);  
                                       bw4=bwareaopen(bw3,3)+oldbw;
                                       oldbw=bw4;
                                   end 
                                   oldbw(oldbw>1)=1;
                                   cc=bwconncomp(oldbw);
                                    Centroid = regionprops3(cc, 'Centroid','Volume',"PrincipalAxisLength");
                                    temptotal=[temptotal size(Centroid.Centroid,1)];
                                    idx=ones(1,size(Centroid.Centroid,1));
                                     high=max(data(:));
                                        low=min(data(:));
%                                     figure;
                                    sn=[];
                                    signal=[];
                                    if ~isempty(Centroid.Centroid)
                                    for i=1:size(Centroid.Centroid,1)
                                        a1=int16(Centroid.Centroid(i,1));
                                        b1=int16(Centroid.Centroid(i,2));
                                        c1=int16(Centroid.Centroid(i,3));

                                        if a1-5<=0 | a1+5>x2-x1 | b1-5<=0 | b1+5>y2-y1 | a1-1<=0 | a1+1>x2-x1 | b1-1<=0 | b1+1>y2-y1
                                            sn(i)=0;
                                        else
                                            s=max(max(data(b1-1:b1+1,a1-1:a1+1,c1)));
                                            signal(i)=s;
                                        n=min([data(b1-5,a1,c1),data(b1+5,a1,c1),data(b1,a1+5,c1),data(b1,a1-5,c1)]);
                                        sn(i)=s/n;
                                        end
                                    end
%                                      disp(sum(idx));
                                   
                                   bitinfo(l,bit).x=Centroid.Centroid(:,1);
                                   bitinfo(l,bit).y=Centroid.Centroid(:,2);
                                   bitinfo(l,bit).z=Centroid.Centroid(:,3)+shift(3); 
                                   bitinfo(l,bit).s=signal;
                                   bitinfo(l,bit).snratio=sn;
                                   tempcount=[tempcount sn];  
%                                    hold on;
%                                    scatter(bitinfo(l,bit).x+x1,bitinfo(l,bit).y+y1,'ro');
%                                    hold on;
                                   
                   end
                end
            end   
            disp(['finished spotfinding for Bit ',num2str(bit), ' average preselected spots per cell:' num2str(mean(temptotal)) ' average signal-noise ratio:' num2str(mean(tempcount)) 'for FOV:', num2str(fovID) ]);
            toc;
        end           
             parsave_ChIPFISH_signal(strcat(mapPath,['FOV_' num2str(fovID),'data.mat']),cellBoundaryStruct,ValidCell,bitinfo,shiftall); 
    end 
    end              
           close all;
    end    
end

%%
marker='o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*o^><+s*';
histFcn = @(x) histogram(x, 'Normalization', 'pdf');
for fovID=1:100
%% find the spot and barcode 
   if ~exist(strcat(mapPath,['spots_FOV' num2str(fovID) '.png']),'file') && exist(strcat(mapPath,['FOV_' num2str(fovID),'data.mat']),'file')  
     threshold=1.4;
      temp=load([mapPath,'FOV_' num2str(fovID),'data.mat']);
      cellBoundaryStruct=temp.cellBoundaryStruct;
      bitinfo=temp.bitinfo;
      ValidCell=temp.ValidCell;
      shiftall=temp.shiftall;
      spotlocalization={};
       spotnumberpercell=[];
       decodeinfo=[];
                 errorinfo=[];
                 errorbit={};
       spoton=zeros(length(cellBoundaryStruct),200,bitnumber);
       for l=1:size(bitinfo,1) %length(cellBoundaryStruct)
           if ValidCell(l)
               allspot=[];
               spotbit=[];
                for bit=1:size(bitinfo,2)%bitnumber
                    idx=bitinfo(l,bit).snratio'>=threshold;
                    allspot=[allspot;[bitinfo(l,bit).x(idx),bitinfo(l,bit).y(idx),bitinfo(l,bit).z(idx)*2]];
                    spotbit=[spotbit, ones(1,length(bitinfo(l,bit).x(idx)))*bit];

%                     end
                end
                      idx=dbscan(allspot,2.5,3);
                           c=unique(idx);
                          
                           id=idx==-1;
                    disp(['Expected decoded spot number for cell ' num2str(l),':' num2str(length(allspot)/4),' Actual Cluster number:', num2str(length(c)-1)]);
                 
                 center=[];
                 index=0;
                 for i=1:length(c)-1
                     id=idx==i;
                   %if inpolygon(x1+mean(allspot{l}(id,1)),y1+mean(allspot{l}(id,2)),cellBoundaryStruct(l).cytoplasmBoundary(:,2),cellBoundaryStruct(l).cytoplasmBoundary(:,1))==1 && std(allspot{l}(id,1))>2
                   barcode=zeros(1,bitnumber);
                       barcode(spotbit(id))=1;
                   if sum(id)<=5
                       dist2=pdist2(code,barcode);
                       err=min(dist2);
                       idx1=find(dist2==err);
                       if err<=1
                          index=index+1; 
                          center(index,1)=mean(allspot(id,1));
                          center(index,2)=mean(allspot(id,2));
                          center(index,3)=mean(allspot(id,3));
                          decodeinfo(l,index)=idx1;
                          errorinfo(l,index)=err;
                          errorbit{l,index}=find(abs(code(idx1,:)-barcode));
                       end
                   else
                       
                       idx3=kmeans(allspot(id,:),ceil(sum(id)/5));
                       idx3=dbscan(allspot(id,:),1.5,3);
%                        figure;
                       color='rgbky';
                       idnum=find(id);
                       w=0;
                       for k=1:length(unique(idx3))-1%ceil(sum(id)/5)
%                            scatter3(allspot(idnum(idx3==k),1),allspot(idnum(idx3==k),2),allspot(idnum(idx3==k),3),[],color(k));
%                            hold on;
                           barcode=zeros(1,bitnumber);
                           barcode(spotbit(idnum(idx3==k)))=1;
                           dist2=pdist2(code,barcode);
                               err=min(dist2);
                               idx1=find(dist2==err);
                               if err<=1
                                  w=w+1; 
                                  index=index+1; 
                                  center(index,1)=mean(allspot(idnum(idx3==k),1));
                                  center(index,2)=mean(allspot(idnum(idx3==k),2));
                                  center(index,3)=mean(allspot(idnum(idx3==k),3));
                                  decodeinfo(l,index)=idx1;
                                  errorinfo(l,index)=err;
                                  errorbit{l,index}=find(abs(code(idx1,:)-barcode));
                               end
                       end
                      % disp(['find overlapping spots,Found ' num2str(w) ' out of ' num2str(ceil(sum(id)/5))]);
                      disp(['find overlapping spots,Found ' num2str(w) ' out of ' num2str(sum(id)) 'spots, out of ' num2str(length(unique(idx3))-1) 'clusters.' ]);
                   end
                 end 
                 spotnumberpercell(l)=index;    
                 spotlocalization{l}=center;
           end
       end

%%
   disp(['Start plot the spots for FOV',num2str(fovID)]);
figHandle = figure(...
            'Name', ['spots_FOV' num2str(fovID)], ...
            'visible', 'on');
       for i=1:length(cellBoundaryStruct)
        plot(cellBoundaryStruct(i).cytoplasmBoundary(:,2),cellBoundaryStruct(i).cytoplasmBoundary(:,1),'k');
        hold on;
         end
             
            hold on;
            c='bgrcmyk';
            s='ox+*sdv^<>ph';
            for l=1:size(bitinfo,1)
                if ValidCell(l)& ~isempty(spotlocalization{l})
                    x1=max(min(cellBoundaryStruct(l).cytoplasmBoundary(:,2))-3,1);
                       y1=max(min(cellBoundaryStruct(l).cytoplasmBoundary(:,1))-3,1);
                       temp=spotlocalization{l};
                       temp1=errorinfo(l,:); 
                       temp2=decodeinfo(l,:);
                      
                    for i=1:min(spotnumberpercell(l),length(temp1))
                        label='';
                        if temp1(i)<=1 & temp2(i)>0
                               tf='none';
                               if temp1(i)==0 && temp2(i)>0
                                  tf='filled';
                               end
                               if temp2(i)>127 % number of target loci
                                   label='k.';
                               elseif temp2(i)>0
                                   chr=str2num(Name{temp2(i)}(4:end));
                                   label=[c(1+mod(chr,7)),s(1+mod(chr,12))];

%                                else 
%                                    label='ko';
                               end
                               if strcmp(tf,'none')
                               scatter(temp(i,1)+x1,temp(i,2)+y1,label);
                               else
                               scatter(temp(i,1)+x1,temp(i,2)+y1,label,tf);
                               end
                               hold on;
                        end
                    end
                end
            end
         SaveFigure(figHandle, 'overwrite', true, ...
                'formats', {'fig', 'png'}, ...
                'savePath', mapPath);               
             parsave_ChIPFISH(strcat(mapPath,['FOV_' num2str(fovID),'data.mat']),spotlocalization,spotnumberpercell,spoton,decodeinfo,errorinfo,errorbit,cellBoundaryStruct,ValidCell,bitinfo,shiftall); 
        %%
                    
           close all;
   end
end
%% load data

total=0;
fov=0;
localization={};
numberpercell=[];
barcode=[];
decode={};
error={};
errorbit={};
cellnum=0;
cellID=[];
FOV=[];
badfov=[];
cellmatrix=zeros(4000,366);
for fovID=1:100
    if exist(strcat(mapPath,['FOV_' num2str(fovID),'data.mat']), 'file') == 2 && sum(badfov==fovID)==0
        temp1=load(strcat(mapPath,['FOV_' num2str(fovID),'data.mat']));
%         temp2=load(strcat(mapPath1,['FOV_' num2str(fovID),'data.mat']));
        idx=boolean(temp1.ValidCell);
        idx=idx(1:size(temp1.decodeinfo,1));
        fov=[fov,fovID*ones(1,sum(idx))];
        localization{fovID}=temp1.spotlocalization(idx);
        numberpercell{fovID}=temp1.spotnumberpercell(idx);
        bound=temp1.cellBoundaryStruct(idx);
        decode{fovID}=temp1.decodeinfo(idx,:);
        error{fovID}=temp1.errorinfo(idx,:);       
         errorbit{fovID}=temp1.errorbit(idx,:);   
        for i=1:sum(idx)            
            idx1=zeros(1,numberpercell{fovID}(i));
            temp2=error{fovID}(i,1:numberpercell{fovID}(i));           
                    cellnum=cellnum+1;
                    FOV(cellnum)=fovID;
                    cellID(cellnum)=i;
                    cellbound{cellnum}=bound(i).cytoplasmBoundary;
                    temp3=decode{fovID}(i,1:numberpercell{fovID}(i));
                    for j=1:366
                            cellmatrix(cellnum,j)=sum(temp3==j);
                    end
        end    
    end
end

%%
temp=sum(cellmatrix,2);
fig = figure(...
            'Name', ['count per cell'], ...
            'visible', 'on');
histogram(temp);        
title(['Median decoded number per cell: ' num2str(median(temp))]);
SaveFigure(fig, 'overwrite', true, ...
                'formats', {'fig', 'png'}, ...
                'savePath', decodePath); 
%%
clusternum=length((unique(c)-1));           
color=jet(clusternum);
fig = figure(...
            'Name', ['UMAP'], ...
            'visible', 'on');
for i=1:clusternum
    scatter(a(c==i,1),a(c==i,2),[],color(i,:),'.');
    hold on;
end
xlabel('UMAP-X');
ylabel('UMAP-Y');
SaveFigure(fig, 'overwrite', true, ...
                'formats', {'fig', 'png'}, ...
                'savePath', decodePath); 
            disp(clusternum)
%%

%% measure the number of spots exact match to the codebook
numberpercellzero=[];
numberpercellone=[];
errorzero=[];
errorone=[];
errortwo=[];
biterror=[];
for fovID=1:100      
   if exist(strcat(mapPath,['FOV_' num2str(fovID),'data.mat']), 'file') == 2 
    for i=1:size(error{fovID},1)
        temp=error{fovID}(i,1:numberpercell{fovID}(i));
        idx0=find(temp==0);
        errorzero=[errorzero,decode{fovID}(i,idx0)];     
        idx1=find(temp==1);
         errorone=[errorone,decode{fovID}(i,idx1)];
         biterror=[biterror,errorbit{fovID}{i,idx1}];
         
    end       
   end
end
fig = figure(...
            'Name', ['error ratio'], ...
            'visible', 'on');
    bar([length(errorzero),length(errorone),length(errortwo)]);
    title(['Exact ratio: ' num2str(length(errorzero)/(length(errorzero)+length(errorone)))]);
    xlabel('Error number');
    ylabel('Loci number');
    set(gca, 'XTickLabel', {0,1,2})
SaveFigure(fig, 'overwrite', true, ...
                'formats', {'fig', 'png'}, ...
                'savePath', decodePath);  
            
fig = figure(...
            'Name', ['bit error'], ...
            'visible', 'on');   
e=[];
for bit=1:bitnumber
    e(bit)=sum(biterror==bit);
end
bar(e);
SaveFigure(fig, 'overwrite', true, ...
                'formats', {'fig', 'png'}, ...
                'savePath', decodePath); 
%% plot target loci vs blank

all=[errorone,errorzero];
decodeall=[];
for i=1:366
    decodeall(i)=sum(all==i);
end

fig = figure(...
            'Name', ['barcode_abundance'], ...
            'visible', 'on');  
[sortednA, sind] = sort(decodeall/length(fov), 'descend');
    x = decodeall(sind);
    localBlankInds = find(sind>127);
    bar(1:length(sortednA), sortednA, 1, 'b', 'EdgeColor', 'none'); hold on;
    bar(localBlankInds, sortednA(localBlankInds), 1, 'r', 'EdgeColor', 'none');
    set(gca,'YScale', 'log');
    xlabel('Barcode ID');
    ylabel('Counts');
%     xlim([0 141]);    
    maxBlank = max(decodeall(127:end));
    title(['Number above: ' num2str(sum(decodeall(1:127) > maxBlank))]);
  SaveFigure(fig, 'overwrite', true, ...
                'formats', {'fig', 'png'}, ...
                'savePath', decodePath);   
save([decodePath 'decode'],'decodeall','fov','localization','numberpercell','barcode','decode','error','cellmatrix','cellID','FOV','cellnum','cellbound');  
%%
pos= dlmread('\\rembrandt\Data\ChIP_FISH\setting\Position_100_new_60x.txt');
all=zeros(int64((max(pos(:,1))-min(pos(:,1)))/200*2048),int64((max(pos(:,2))-min(pos(:,2)))/200*2048));
for fovID=1:length(pos)
    a=ReadDax([DAPIPath,'zscan_',num2str(fovID-1,'%03d'),'_0.dax'], 'startFrame', 5, 'endFrame', 5,'verbose',false);
    all(int16((pos(fovID,1)-min(pos(:,1)))/200*2048+1):int16((pos(fovID,1)-min(pos(:,1)))/200*2048+2048),int16((pos(fovID,2)-min(pos(:,2)))/200*2048)+1:int16((pos(fovID,2)-min(pos(:,2)))/200*2048+2048))=flip(int16(double(a)./ratio405/10),1);    
end
figure;
imshow(all,[1.2*min(all(:)),0.8*max(all(:))]);
%%
%% plot spots
for b=1:140
    figure;
    imshow(all,[1.4*min(all(:)),0.6*max(all(:))]);
    hold on;
    for fovID=1:length(pos)
        if exist(strcat(mapPath,['FOV_' num2str(fovID),'data.mat']),'file') 
            temp=load(strcat(mapPath,['FOV_' num2str(fovID),'data.mat']));
            for l=1:size(temp.errorinfo,1)
                       x1=max(min(temp.cellBoundaryStruct(l).cytoplasmBoundary(:,2))-3,1);
                       y1=max(min(temp.cellBoundaryStruct(l).cytoplasmBoundary(:,1))-3,1);
                for j=1:temp.spotnumberpercell(l)
                     if temp.decodeinfo(l,j)==b 
                         scatter(temp.spotlocalization{l}(j,1)+x1+int16((pos(fovID,2)-min(pos(:,2)))/200*2048)+1,2048-temp.spotlocalization{l}(j,2)-y1+int16((pos(fovID,1)-min(pos(:,1)))/200*2048+1), 25,'ro','filled');
                         hold on;
                     end   
                end 
            end
        end
    end
    title(Gene{b});
     saveas(gcf,[decodePath 'Loci_' num2str(b) '_' Gene{b} '.png']);
     close all;
end


%% plot positive cells in the cortex
%Candidate={'Satb2','Cux2','Nr2f1','Rorb','Fezf2','Bcl11b','Tbr1','Sox5','Erbb4','Olig2','Sst','Pou3f2','Slc17a6','Dlx1','Reln','Lamp5'};
Candidate={'Gpr88'};
plotPath=SetFigureSavePath([analysisSavePath,'\individual_plot\'], ...
     'makeDir', true);
c='rbyg';
b=[];
badfov=[6,7,8,9,10,11]; % fovs that don't have any cells for the tissue slices
for gene=1:length(Candidate)
    figure;
        for i=1:cellnum
            if sum(badfov==FOV(i))==0
                fill(cellbound{i}(:,2)+(pos(FOV(i),2)-min(pos(:,2)))/200*2048+1,2048-cellbound{i}(:,1)+(pos(FOV(i),1)-min(pos(:,1)))/200*2048+1,[0.8,0.8,0.8],'LineStyle','none');
            hold on;
            end
       end
       b=find(contains(Gene,Candidate{gene}));

            for fovID=1:length(pos)
                if exist(strcat(mapPath,['FOV_' num2str(fovID),'data.mat']),'file') && sum(badfov==fovID)==0
                    temp=load(strcat(mapPath,['FOV_' num2str(fovID),'data.mat']));
                    for l=1:size(temp.errorinfo,1)
                               x1=max(min(temp.cellBoundaryStruct(l).cytoplasmBoundary(:,2))-3,1);
                               y1=max(min(temp.cellBoundaryStruct(l).cytoplasmBoundary(:,1))-3,1);
                        for j=1:temp.spotnumberpercell(l)
                             if temp.decodeinfo(l,j)==b(1)
                                 scatter(temp.spotlocalization{l}(j,1)+x1+int16((pos(fovID,2)-min(pos(:,2)))/200*2048)+1,2048-temp.spotlocalization{l}(j,2)-y1+int16((pos(fovID,1)-min(pos(:,1)))/200*2048+1), 10,'ro','filled');
                                 hold on;
                             end   
                        end 
                    end
                end
            end
     axis off
     camroll(180);
     title(Candidate{gene});
     saveas(gcf,[plotPath Candidate{gene} '_thre1.4.png']);
end
close all;


