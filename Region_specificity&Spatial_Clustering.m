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

%% define fov
cortex=[423,394,393,392,391,383,384,349,350,348,342,343,344,347,348,349,385,386,389,390,424,425,426,427,428,431,432,433,467,468,469,470,474,475,476,477,508,509,510,517,518,519,520,549,550,551,548,521,522,547,563,564,565,566];
Striatum=[398,399,422,421,420,419,418,417,416,443,442,441,440,439,438,437,436,465,464,463,462,461,460,459,458,485,484,483,482,481,480,479,478,506,505,504,503,502,501,500,527,526,525,524,523,522,546,545,544,543,542,568,567];
Diencephalon=[395,396,397,358,357,356,355,354,353,352,351,350,341,340,339,338,337,336,335,334,333,316,315,314,313,312,311,310,309,308,307,298,297,296,295,294,293,292,291,290,275,274,273,272,271,270,269,268,267,446,445,444,415,414,413,404,403,402,401,400,376,377,378,379,380,381,375,374,373,372,371,362,361,360,359,358,333,332,331,330,329,319,318,317,290,289,288];
forebrain=[cortex,Striatum,Diencephalon];
hindbrain=[62,63,64,65,66,67,68,69,70,71,32,33,34,35,36,74,75,76,77,78,79,80,81,82,83,105,106,107,108,109,110,111,112,113,116,117,118,119,120,121,122,123,124,146,147,148,149,150,151,152,153,154,155,158,159,160,161,162,166,167,168,169,170,185,186,187,188,189,190,191,193,194,195,196,197,200,201,202,203,204,205,206,207,208,209,210,211,212,227,228,229,230,231,232,233,234,235,236,237,238,243,244,245,246,247,248,249];
midbrain=[266,300,299,266,265,264,263,260,259,258,257,256,255,254,227,226,225,224,223,222,221,218,217,216,215,214,213,184,183,182,181,180,179,176,175,174,173,172,171,142,141,140,139,138,137,134,133,132,131,130,129,128,127,102,103,104,101,100,99,98,97,89,88,87,86,85,84,61,60,125];
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
cgo = clustergram(regioncount','ColumnPDist','Correlation','Linkage','ward','Colormap','redbluecmap','Standardize','Column');
set(cgo,'Dendrogram',10);
%%
Y=pdist(regioncount);
z=linkage(Y);
figure;
dendrogram(z)
