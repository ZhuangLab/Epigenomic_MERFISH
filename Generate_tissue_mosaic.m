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
shiftPath=SetFigureSavePath([analysisSavePath,'\shift\'], ...
     'makeDir', true);
pos= dlmread('');
%% plot decoded spots with colors corresponding to their identity
color=hsv(127);
idx=randperm(127);
color=color(idx,:);
figure;
for i=1:127
    scatter(i,0,300,color(i,:),'s','filled');
    hold on;
end
%% define fov
cortex=[11,12,13,14,23,26,25,27,22,21,78,77,65,66,64,63,84,87,86,89,97,96,95,115,147,146,116,119,120,335,338];
caudal=[76,67,75,68,74,69,61,62,88,90,91,73,28,29,32,70,60,92,93,94,148,152,153,59,71,72,47,48,58,154,151,149,144,143,150,156,155,57,49,50,56,168,167,157,142,141,158,166,169,170,165,159];
thylamus=[7,8,5,6,342,343,346,1,2,34,31,30,347,356,370,371,35,340,341,344,369,368,38,37,39,40,45,50,345,348,36];
hypothylamus=[182,51,52,53,42,43,44,41,180,179,178,181,177,187,186,183,366,365,361,364,185,191,359,358,354,355,362];
hindbrain=[211,212,213,218,217,214,215,216,219,220,248,251,252,249,285,288,289,290,246,245,243,240,225,224,227,226,239,236,228,229,235,231,232,237,238,241,234,286,287,244,351];
midbrain=[255,252,253,254,258,259,260,262,261,268,267,270,280,281,266,265,250,284,283,282,316,317,318,313,314,315,294,293,295,296,297,312,298,299,300,309,305,303,304,303,302,340,307,308,309,310,306,303];
%% plot DAPI
all=zeros(int64((max(pos(:,1))-min(pos(:,1)))/210*2048),int64((max(pos(:,2))-min(pos(:,2)))/210*2048));
for fovID=1:length(pos)
    a=ReadDax([DAPIPath,'Epi-750s1-650s1-560s1-488s1-405s1_',num2str(fovID-1,'%03d'),'_0.dax'], 'startFrame', 5, 'endFrame', 5,'verbose',false);
    all(int64((pos(fovID,1)-min(pos(:,1)))/210*2048+1):int64((pos(fovID,1)-min(pos(:,1)))/210*2048+2048),int64((pos(fovID,2)-min(pos(:,2)))/210*2048)+1:int64((pos(fovID,2)-min(pos(:,2)))/210*2048+2048))=flip(a,2);    
    
end
figure;
imshow(all,[1.2*min(all(:)),0.1*max(all(:))]);
hold on;
for fovID=1:length(pos)
text(double((pos(fovID,2)-min(pos(:,2)))/210*2048+1024),double((pos(fovID,1)-min(pos(:,1)))/210*2048+1024),num2str(fovID),'Color','w','FontSize',10);
hold on;
end
%% plot cells
figure;
imshow(all,[1.2*min(all(:)),0.1*max(all(:))]);
hold on;
regioncount=zeros(366,6);
badfov=[276,275,326,327,328,105,106,109,110,323,322,17,18,104,107,108,325,324,321,329,81,103,101,100,111,112,274,277,273,272,113,114,117,122,123,129,131,136,138,139,162,176,188,189,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,221,222];
% badfov is the fov with no cells for the tissue
for fovID=1:length(pos)
    tempcount=zeros(366,1);
        if exist(strcat(mapPath,['FOV_' num2str(fovID),'data.mat']),'file')&& sum(badfov==fovID)==0
            temp1=load(strcat(mapPath,['FOV_' num2str(fovID),'data.mat']));
            if sum([145,146,125,126,133,132]==fovID)>0
                temp1=load(strcat(mapPath,['FOV_' num2str(150),'data.mat']));
            end
            local=temp1.spotlocalization;
            E=temp1.errorinfo;
            D=temp1.decodeinfo;
                for i=1:length(D)
                    if contains(Name{D(i)},['chr'])
                       c=str2num(Name{D(i)}(4:end));
                       scatter(2048-local(i,1)+(pos(fovID,2)-min(pos(:,2)))/210*2048+1,local(i,2)+(pos(fovID,1)-min(pos(:,1)))/210*2048+1,3,color(c,:),'filled');
                       hold on;
                    end
                 end
        end
      disp(fovID)
end


