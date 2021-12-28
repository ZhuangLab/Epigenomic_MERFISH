
load(''); % single cell profile include localization of the cells and the epigenomic profiles in each cell
 
 %%
%spots=[139.68,4579.4715;8692.77,6661.66;18158.80,10434];
center=[-11559.09 , 71243.03];
locidepth=zeros(length(cellbound),2);
lociname=strings(10,length(cellbound));
[~,~,codebook]=xlsread('');
code=cell2mat(codebook(1:366,5:28));
Name=codebook(1:366,4);
%% cell density cross layers
celldis=[];
for i=1:length(cellbound)
   celldis(i)=pdist([[mean(cellbound{i}(:,2)),mean(cellbound{i}(:,1))];center]);
end
celldis=(celldis-min(celldis))/(max(celldis)-min(celldis));

%% color the cells with different radial distance
figure;
for i=1:length(cellbound)
    color='r';
    if celldis(i)>0.82
        color='r';
    elseif celldis(i)>0.65
        color='m';
    elseif celldis(i)>0.45
        color='b';
    else
        color='k';
    end
    plot(cellbound{i}(:,2),cellbound{i}(:,1),color);
    hold on;
end
camroll(180)
%%
Newdepth=zeros(127,4);
pvalue=zeros(1,127);
oncells=zeros(127,4);
offcells=zeros(127,4);
cellcount=zeros(length(ida),4);
allcell=[];
total=[sum(celldis>0.82),sum(celldis>0.65 & celldis<=0.82),sum(celldis>0.45 & celldis<=0.65),sum(celldis<=0.45)];
for i=1:127
    temp=celldis(newcellmatrix(:,i)>0);
    allcell(i)=length(temp);
   x1=[sum(temp>0.82),sum(temp>0.65 & temp<=0.82),sum(temp>0.45 & temp<=0.65),sum(temp<=0.45)];
   x2=[sum(celldis>0.82)-sum(temp>0.82),sum(celldis>0.65 & celldis<=0.82)-sum(temp>0.65 & temp<=0.82),sum(celldis>0.45 & celldis<=0.65)-sum(temp>0.45 & temp<=0.65),sum(celldis<=0.45)-sum(temp<=0.45)];
   %x2=[sum(celldis>0.82),sum(celldis>0.65 & celldis<=0.82),sum(celldis>0.45 & celldis<=0.65),sum(celldis<=0.45)];
   oncells(i,:)=x1;
   offcells(i,:)=x2;
   for j=1:4
       Newdepth(i,j)=x1(j)/total(j);
   end
   pvalue(i)= chi2test([x1',x2']);
end
imaging=Newdepth;
%%
newcount2=imaging;
for i=1:length(imaging)
    newcount2(i,:)=(imaging(i,:)-mean(imaging(i,:)))./std(imaging(i,:));
end
[~,maxid]=max(newcount2,[],2);
index=1;
name=[];
label=[];
pvalue1=[];
sortedregioncount1=newcount2;
for i=1:4
    idx2=maxid==i;
    [temp,idx3]=sortrows(newcount2(idx2,:),[i],'descend');
    sortedregioncount1(index:index+sum(idx2)-1,:)=temp;
    index=index+sum(idx2);
    tempid=find(idx2);
    name=[name;tempid(idx3)];
    pvalue1=[pvalue1,pvalue(tempid(idx3))];
    label=[label;Name(tempid(idx3))];
end

%%
figure;
heatmap(label,{'II/III','IV','V','VI'},sortedregioncount1','Colormap',redbluecmap);
grid off;
figure;
pvalue2=pvalue1;
pvalue2(pvalue2>0.2)=0.25;
heatmap(label,{'P-value'},pvalue2,'Colormap',jet);
title ('all cells');


