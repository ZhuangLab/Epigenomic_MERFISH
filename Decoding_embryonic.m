%% Measure the MERFISH and cell identities
%% Tian LU
%% 3/1/2017
%% Setup path and parameters
% Define data path
dataPath='';
analysisSavePath = SetFigureSavePath([dataPath,'\FOV\'], ...
     'makeDir', true);
RawPath = ''; % Rawdata
DAPIPath=''; % dapi images for segmentation
% Useful data structure for spotfinding
all=readtable(''); % readout path

[~,~,codebook]=xlsread(''); %Codebook path
code=cell2mat(codebook(1:366,4:27));
Name=codebook(1:366,1);
color=all.color;
desiredround=all.imagingRound;
bitnumber=24;
roundnumber=12;
znumber=30;

mapPath=SetFigureSavePath([analysisSavePath,'\loci_localization\'], ...
     'makeDir', true);
decodePath=SetFigureSavePath([analysisSavePath,'\decode_per_cell\'], ...
     'makeDir', true);
%% ------------------------------------------------------------------------
% Start logging
%%-------------------------------------------------------------------------

%% Correct the illumination of 647, 750, 561,405

if exist([dataPath 'illumination.mat'],'file')
    load([dataPath 'illumination']);
else
    Path=[RawPath];
    file=['Epi-750s30-650s30-555s30-488s30_(?<fov>[0-9]+)_0'];    
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
    
  %  
    Path=DAPIPath;
    file=['Epi-750s1-650s1-560s1-488s1-405s1_(?<fov>[0-9]+)_0'];    
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

   save([dataPath 'illumination'],'ratio750','ratio647','ratio405');

end  

%% Create parallel pool
if isempty(gcp('nocreate'))
    p = parpool(20); % Set this number to control the number of parallel workers.
                     % Polite usage would suggest these maximum values: morgan = 20, cajal = 10
else
    p = gcp;
end


%%
pos= dlmread(''); % list of X Y position for each fov

histFcn = @(x) histogram(x, 'Normalization', 'pdf');
tform1pre=[];
tform2pre=[];
load([dataPath 'illumination']);
badfov=[];
parfor fovID =1:length(pos)
    if ~exist(strcat(mapPath,['FOV_' num2str(fovID),'data.mat']),'file')
       %% load raw images
       warning('off','all')
        % start finding spots in 3D
        disp(['Start mapping the spots to cells: FOV',num2str(fovID)]);       
        bitinfo=[];       
        shiftall={};
       %map spot based on 561 beads 
        fixed=ReadDax([RawPath,'Epi-750s30-650s30-488s30_',num2str(fovID-1,'%03d'),'_01.dax'], 'startFrame', znumber*2+1, 'endFrame', znumber*3,'verbose',false);
        data=fixed;
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
       shiftpre=[0,0,0];
       %% load raw images
       warning('off','all')
        % start finding spots in 3D
        disp(['Start mapping the spots to cells: FOV',num2str(fovID)]);       
        bitinfo=[];       
        for bit=1:bitnumber
                        tic;
            shift=[0,0,0];          
            switch color(bit)
                case 750
                    a=ReadDax([RawPath,'Epi-750s30-650s30-488s30_',num2str(fovID-1,'%03d'),'_', num2str(desiredround(bit),'%02d'),'.dax'], 'startFrame', 1, 'endFrame', znumber,'verbose',false);                  
                    a=double(a)./ratio750/100;
                case 650
                    a=ReadDax([RawPath,'Epi-750s30-650s30-488s30_',num2str(fovID-1,'%03d'),'_', num2str(desiredround(bit),'%02d'),'.dax'], 'startFrame', znumber+1, 'endFrame', znumber*2,'verbose',false);  
                    a=double(a)./ratio647/100;
                    a=flip(a,3);
            end
            if bit>1 && desiredround(bit)==desiredround(bit-1) 
                   shift=shiftpre;
%                     disp('fov:');
%                     disp(fovID);
            elseif bit>1
              moving=ReadDax([RawPath,'Epi-750s30-650s30-488s30_',num2str(fovID-1,'%03d'),'_', num2str(desiredround(bit),'%02d'),'.dax'], 'startFrame', znumber*2+1,'endFrame',znumber*3,'verbose',false);   

                 data=moving;
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
            tic;
             x=[];
             y=[];
             z=[];
                     for focal=1:30
            [data,I]=max(a(:,:,focal),[],3);
%             figure;
            h=histogram(data);
            ccum=zeros(1,h.NumBins);
            for i=1:h.NumBins
            ccum(i)=sum(h.BinCounts(1:i));
            end            

            idx=find(ccum>0.99999*ccum(end));  
           bw=data>h.BinEdges(idx(1));   
           %%bw=data>0.9*max(data(:));     
           bw2=imfill(bw,'holes');
           bw3=bw2-bwareaopen(bw2,100); 
           bw4=bwareaopen(bw3,3);
           oldbw=bw4;
           thre=0.99;
           for threshold=0.99995:-0.0001:thre%bitthreshold(bit)
               idx=find(ccum>threshold*ccum(end));  
           bw=data>h.BinEdges(idx(1));
               %bw=data>threshold*max(data(:));     
               bw2=imfill(bw,'holes');
               bw3=bw2-bwareaopen(bw2,100);  
               bw4=bwareaopen(bw3,5)+oldbw;
               oldbw=bw4;
           end 
           oldbw(oldbw>1)=1;
           cc=bwconncomp(oldbw,4);
%            figure;
%            imshow(oldbw);
            Centroid = regionprops(cc, 'Centroid');
            temp1=[];
            temp2=[];
            idx=zeros(1,length(Centroid));
%             figure;
%            imshow(imadjust(int16(data)));
%            hold on;
            for i=1:length(Centroid)
                a1=int16(Centroid(i).Centroid(1));
                b1=int16(Centroid(i).Centroid(2));
%                 scatter(a1, b1,50,'ro');
                if a1-5<=0 | a1+5>2048 | b1-5<=0 | b1+5>2048 | a1-1<=0 | a1+1>2048 | b1-1<=0 | b1+1>2048
                                            sn(i)=0;
                                        else
                                            s=max(max(data(b1-1:b1+1,a1-1:a1+1)));
                                        n=min([data(b1-5,a1),data(b1+5,a1),data(b1,a1+5),data(b1,a1-5)]);
                                        sn=s/n;
                end
                                        if sn>1.4
                                            idx(i)=1;
                                        end
                temp1(i)=a1;
                temp2(i)=b1;
            end
          idx=boolean(idx);
           x=[x,temp1(idx)];
          y=[y,temp2(idx)];
          z=[z,focal*ones(1,sum(idx))]; 
           
                     end
           toc;
           %%
           idx=dbscan([x',y'],0.5,2);
           fx=[];
           fy=[];
           fz=[];
           for i=1:length(unique(idx))-1
               id=idx==i;
               fx(i)=mean(x(id))+shift(1);
               fy(i)=mean(y(id))+shift(2);
               fz(i)=mean(z(id))+shift(3);
           end
           bitinfo(bit).x=fx;
           bitinfo(bit).y=fy;
           bitinfo(bit).z=fz;
           disp(['finished spotfinding for Bit ',num2str(bit), ' average preselected spots per cell:' num2str(length(x)) ' average spot per cell:' num2str(length(bitinfo(bit).x)) ' for FOV ',num2str(fovID)]); 
           toc;
        end
    
%% find the spot and barcode 
        decodeinfo=[];
                 errorinfo=[];
                 errorbit={};
       spoton=zeros(200,bitnumber);    
      spotlocalization={};
       spotnumberpercell=[];
       spoton=zeros(10000,bitnumber);
               allspot=[];
               spotbit=[];
                for bit=1:bitnumber
                    allspot=[allspot;[bitinfo(bit).x',bitinfo(bit).y',bitinfo(bit).z']];
                    spotbit=[spotbit, ones(1,length(bitinfo(bit).x))*bit];
                end
                      idx=dbscan(allspot,2.5,3);
                           c=unique(idx);
                           center=zeros(length(c)-1,2);
                           id=idx==-1;
                disp(['Expected decoded spot number for cell :' num2str(length(allspot)/4),' Actual Cluster number:', num2str(length(c)-1)]);
                 
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
                          decodeinfo(index)=idx1;
                          errorinfo(index)=err;
                          errorbit{index}=find(abs(code(idx1,:)-barcode));
                       end
                   else
                       
                       idx3=dbscan(allspot(id,:),1.5,3);
%                        figure;
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
                                  decodeinfo(index)=idx1;
                                  errorinfo(index)=err;
                                  errorbit{index}=find(abs(code(idx1,:)-barcode));
                               end
                       end
                      disp(['find overlapping spots,Found ' num2str(w) ' out of ' num2str(sum(id)) 'spots, out of ' num2str(length(unique(idx3))-1) 'clusters.' ]);
                   end
                 end 
                 spotnumberpercell=index;    
                 spotlocalization=center;


%%
   disp(['Start plot the spots for FOV',num2str(fovID)]);
figHandle = figure(...
            'Name', ['spots_FOV' num2str(fovID)], ...
            'visible', 'on');
        c='bgrcmykw';
        s='ox+*sdv^<>ph';
                   temp=spotlocalization;
                   temp1=errorinfo(:); 
                   temp2=decodeinfo(:);
                for i=1:min(spotnumberpercell,length(temp1))
                    if temp1(i)<=1 && temp2(i)>0
                       tf='none';
                       if temp1(i)==0 && temp2(i)>0
                          tf='filled';
                       end
                       if temp2(i)>147 % the number of target loci
                           label='k.';
                       elseif temp2(i)>0
                           chr=str2num(Name{temp2(i)}(4:end));
                           label=[c(1+mod(chr,8)),s(1+mod(chr,12))];
                       else 
                           label='k.';
                       end
                       if strcmp(tf,'none') 
                       scatter(temp(i,1),temp(i,2),label);
                       else
                       scatter(temp(i,1),temp(i,2),label,tf);
                       end
                       hold on;
                    end
                end
         SaveFigure(figHandle, 'overwrite', true, ...
                'formats', {'fig', 'png'}, ...
                'savePath', mapPath);               
             parsave_ChIPFISH_tissue(strcat(mapPath,['FOV_' num2str(fovID),'data.mat']),spotlocalization,spotnumberpercell,spoton,decodeinfo,errorinfo,errorbit,bitinfo,shiftall);  
    %%
                    
          close all;
    end
end

%% load dapi data and generate mosaic for the whole tissue
%%
all=zeros(int64((max(pos(:,1))-min(pos(:,1)))/210*2048),int64((max(pos(:,2))-min(pos(:,2)))/210*2048));
for fovID=1:length(pos)
    a=ReadDax([RawPath,'Epi-750s30-650s30-488s30_',num2str(fovID-1,'%03d'),'_01.dax'], 'startFrame', 15, 'endFrame', 15,'verbose',false);
    %a=ReadDax([DAPIPath,'Epi-750s1-650s1-560s1-488s1-405s1_',num2str(fovID-1,'%03d'),'_0.dax'], 'startFrame', 5, 'endFrame', 5,'verbose',false);
    all(int64((pos(fovID,1)-min(pos(:,1)))/210*2048+1):int64((pos(fovID,1)-min(pos(:,1)))/210*2048+2048),int64((pos(fovID,2)-min(pos(:,2)))/210*2048)+1:int64((pos(fovID,2)-min(pos(:,2)))/210*2048+2048))=flip(a,2);    
    
end
figure;
imshow(all,[1.2*min(all(:)),0.2*max(all(:))]);
hold on;
for fovID=1:length(pos)
text(double((pos(fovID,2)-min(pos(:,2)))/210*2048+1024),double((pos(fovID,1)-min(pos(:,1)))/210*2048+1024),num2str(fovID),'Color','w','FontSize',10);
hold on;
end

%% plot decoded spots for each gene
plotPath=SetFigureSavePath([analysisSavePath,'\individual_loci2\'], ...
     'makeDir', true);
[~,~,temp]=xlsread(''); % Codebook path
region=temp(2:end,8);
cat=temp(2:end,9);
badfov=[141,142,143,140,289,134,78,73,74];
for b=1:150
    figure;
    imshow(all,[1.4*min(all(:)),0.5*max(all(:))]);
    hold on;
    for fovID=1:length(pos)
        if exist(strcat(mapPath,['FOV_' num2str(fovID),'data.mat']),'file')&& sum(badfov==fovID)==0
        temp=load(strcat(mapPath,['FOV_' num2str(fovID),'data.mat']));
        for j=1:length(temp.errorinfo)
             if temp.decodeinfo(j)==b
                 scatter(2048-temp.spotlocalization(j,1)+int64((pos(fovID,2)-min(pos(:,2)))/210*2048)+1,temp.spotlocalization(j,2)+int64((pos(fovID,1)-min(pos(:,1)))/210*2048+1), 10,'ro','filled');
                 hold on;
             end   
        end    
        end
    end
    title(['Loci ' num2str(b) ' ' region{b} ' ' cat{b}]);
     saveas(gcf,[plotPath 'Loci_' num2str(b) '_' region{b} '.png']);
     close all;
end

%% calculate the total decoded spots for each loci
count=zeros(1,366);
    for fovID=1:100%length(pos)
        temp=load(strcat(mapPath,['FOV_' num2str(fovID),'data.mat']));
        for j=1:length(temp.errorinfo)
            for b=1:366
             if temp.decodeinfo(j)==b
                 count(b)=count(b)+1;
             end  
            end
        end         
    end

%% plot target loci vs blank
fig = figure(...
            'Name', ['barcode_abundance'], ...
            'visible', 'on');  
[sortednA, sind] = sort(count, 'descend');
    x = count(sind);
    localBlankInds = find(sind>127);
    bar(1:length(sortednA), sortednA, 1, 'b', 'EdgeColor', 'none'); hold on;
    bar(localBlankInds, sortednA(localBlankInds), 1, 'r', 'EdgeColor', 'none');
    set(gca,'YScale', 'log');
    xlabel('Barcode ID');
    ylabel('Counts');
%     xlim([0 141]);    
    maxBlank = max(count(127:end));
    title(['Number above: ' num2str(sum(count(1:127) > maxBlank))]);
  SaveFigure(fig, 'overwrite', true, ...
                'formats', {'fig', 'png'}, ...
                'savePath', decodePath); 

