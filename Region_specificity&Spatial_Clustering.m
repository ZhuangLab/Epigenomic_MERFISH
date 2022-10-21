dataPath='';
analysisSavePath = SetFigureSavePath([dataPath,'\FOV\'], ...
     'makeDir', true);
RawPath = ''; 
DAPIPath=''; 
% Useful data structure for spotfinding
all=readtable('');

[~,~,codebook]=xlsread('');
code=cell2mat(codebook(1:366,4:27));
Name=codebook(1:366,1);

color=all.color;
desiredround=all.imagingRound;
bitnumber=24;
roundnumber=8;
znumber=30;

mapPath=SetFigureSavePath([analysisSavePath,'\loci_localization\'], ...
     'makeDir', true);
decodePath=SetFigureSavePath([analysisSavePath,'\decode_per_cell\'], ...
     'makeDir', true);
shiftPath=SetFigureSavePath([analysisSavePath,'\shift\'], ...
     'makeDir', true);
 load([dataPath 'illumination']);
 pos= dlmread('');

%% define region
cortex=[];
Striatum=[];
Diencephalon=[];
forebrain=[cortex,Striatum,Diencephalon];
hindbrain=[];
midbrain=[];
%%
regioncount=zeros(366,3);
regiondapi=zeros(1,3);
for fovID=1:length(pos)
    tempcount=zeros(366,1);
        if exist(strcat(mapPath,['FOV_' num2str(fovID),'data.mat']),'file')
            temp=load(strcat(mapPath,['FOV_' num2str(fovID),'data.mat']));
            for j=1:length(temp.errorinfo)
                if temp.decodeinfo(j)<=366
                 tempcount(temp.decodeinfo(j))=tempcount(temp.decodeinfo(j))+1;
                end
            end    
        end
        tempDAPI=ReadDax([DAPIPath,'zscan_',num2str(fovID-1,'%03d'),'_0.dax'], 'startFrame', 5, 'endFrame', 5,'verbose',false);
     temparea=sum(sum(tempDAPI>600));
   if sum(forebrain==fovID)==1
        regioncount(:,1)=regioncount(:,1)+tempcount;
        regiondapi(1)=regiondapi(1)+temparea/100/1000000; % density:spots/mm2
    elseif  sum(midbrain==fovID)==1
        regioncount(:,2)=regioncount(:,2)+tempcount;
        regiondapi(2)=regiondapi(2)+temparea/100/1000000;
    elseif  sum(hindbrain==fovID)==1    
        regioncount(:,3)=regioncount(:,3)+tempcount;
        regiondapi(3)=regiondapi(3)+temparea/100/1000000;
    end
end

regioncount1=regioncount(1:147,:)./regiondapi;
%%
totalfov=[forebrain,midbrain,hindbrain];
regioncount=zeros(147,length(totalfov));
for fovID=1:length(totalfov)
    tempcount=zeros(147,1);
        if exist(strcat(mapPath,['FOV_' num2str(totalfov(fovID)),'data.mat']),'file')
            temp=load(strcat(mapPath,['FOV_' num2str(totalfov(fovID)),'data.mat']));
            for j=1:length(temp.errorinfo)
                if temp.decodeinfo(j)<=147
                 tempcount(temp.decodeinfo(j))=tempcount(temp.decodeinfo(j))+1;
                end
            end    
        end
        regioncount(:,fovID)=regioncount(:,fovID)+tempcount;
end
%% cluster the enhancer based on fov
cgo = clustergram(regioncount','Linkage','weighted','Colormap','redbluecmap','Standardize','Column');
set(cgo,'Dendrogram',22.3);

%%
Y=pdist(regioncount);
z=linkage(Y);
figure;
dendrogram(z)
