% load data
plotPath='';
alldata=h5read('counts.h5ad','/X');
gene=h5read('counts.h5ad','/var');
gene=gene.index;
cellID=h5read('counts.h5ad','/obs');
cellID=cellID.index;
%% find the overlap gene
[~,~,codebook]=xlsread('');
Name=codebook(1:366,4);
[selectedgene,ida,idb]=intersect(gene,Name);
%% load sample 6
varNames = {'cell','b1','b2','b3','b4','b5','b6','b7'} ;
varTypes = {'char','char','char','char','char','char','char','char'} ;
delimiter = ',';
dataStartLine = 2;
extraColRule = 'ignore';
file=6;
opts=delimitedTextImportOptions('VariableNames',varNames,...
                                'VariableTypes',varTypes,...
                                'Delimiter',delimiter,...
                                'DataLines', dataStartLine,...
                                'ExtraColumnsRule',extraColRule); 
   all=readtable(['\MOp_MERFISH\segmented_cells_mouse1sample' num2str(file) '.csv'],opts);
%% find the cells within slice 313

x1=-2000;
x2=2000;
y1=-100;
y2=2000;
idx=zeros(1,height(all));
for i=2:height(all)
    if length(strfind(all.b5{i},','))>1
  if max(str2num(all.b5{i}))>x1 && min(str2num(all.b5{i}))<x2 && max(str2num(all.b6{i}))>y1 && min(str2num(all.b6{i}))<y2
   idx(i)=1;
  end
    end
end
all1=all(boolean(idx),:);
%% find the cells are in the slice 313
[overlap,idx1,idx3]=intersect(cellID,all1.cell);
idx2=zeros(1,length(cellID));
idx2(idx1)=1;
idx2=boolean(idx2);

%% plot cells
count=10;
threshold=[0,0.1,0.5,0.6,1.1,1.7,2.8,4.5,7.3,11.8,19.1];
 for i=1:length(ida)
    % plot all the cells
    figure;
   for j=1:size(all1,1)
       if length(strfind(all1.b5{j},','))>1
      fill(str2num(all1.b5{j}),str2num(all1.b6{j}),[0.8,0.8,0.8],'LineStyle','none');
      hold on;
       end
   end
  % find the cells that are in sample 6
    idx=alldata(ida(i),:)>0;
    expression=alldata(ida(i),idx & idx2);
    [~,t]=histcounts(expression,threshold); 
    goodcell=cellID(idx & idx2);
    for j=1:length(goodcell)
         index=find(strcmp(all1.cell,goodcell(j)));
         if expression(j)>=max(threshold)
             index2=11;
         else
             index2=min(find(t>expression(j)));
         end
         scatter(mean(str2num(all1.b5{index})),mean(str2num(all1.b6{index})),10,index2-1,'filled');
         hold on;
    end
 % plot the cell by their gene expression level
 colormap(cool(count));
     camroll(180)
    title(selectedgene{i});
    axis off
    cbh=colorbar;
    cbh.TickLabels=num2cell(threshold);
     saveas(gcf,[plotPath selectedgene{i} '.png']);
%     close all
end

%%
center=[1539.83 , 5812.2];
celldis=[];
 for j=1:size(all1,1)
     if length(strfind(all1.b5{j},','))>1
   celldis(j)=pdist([[mean(str2num(all1.b5{j})),mean(str2num(all1.b6{j}))];center]);
     end
end
celldis=(celldis-min(celldis))/(max(celldis)-min(celldis));

%%
figure;
for i=1:length(cellbound)
    color='r';
    if celldis(i)>0.68
        color='r';
    elseif celldis(i)>0.55
        color='m';
    elseif celldis(i)>0.31
        color='b';
    else
        color='k';
    end
    plot(str2num(all1.b5{i}),str2num(all1.b6{i}),color);
    hold on;
end
camroll(180)
%%
alldetected=[];
Newdepth=zeros(length(ida),4);
cellcount=zeros(length(ida),4);
pvalue=zeros(1,length(ida));
total=[sum(celldis>0.68), sum(celldis>0.55 & celldis<=0.68),sum(celldis>0.31 & celldis<=0.55),sum(celldis>0 & celldis<=0.31)];
for i=1:length(ida)
        idx=alldata(ida(i),:)>3;
        alldetected(i)=sum(idx & idx2);
    expression=alldata(ida(i),idx & idx2);
    goodcell=cellID(idx & idx2);
    temp=[];
    for j=1:length(goodcell)
         index=find(strcmp(all1.cell,goodcell(j)));
         temp(j)=celldis(index);
    end
   x1=[sum(temp>0.68), sum(temp>0.55 & temp<=0.68),sum(temp>0.31 & temp<=0.55),sum(temp>0 & temp<=0.31)];
   cellcount(i,:)=x1;
   for j=1:4
       Newdepth(i,j)=x1(j)/total(j);
   end
   pvalue(i)=chi2test([x1',(total-x1)']);
end
imaging=Newdepth;
%%
figure;
for i=1:length(pvalue)
    if pvalue(i)<0.01
        sz=20;
    elseif pvalue(i)<0.05
        sz=10;
    else
        sz=5;
    end
    scatter(1,i,sz,'b','filled');
    hold on;
end
axis off

%%
Newdepth=zeros(length(ida),4);

total=[sum(celldis>0.68), sum(celldis>0.55 & celldis<=0.68),sum(celldis>0.31 & celldis<=0.55),sum(celldis>0 & celldis<=0.31)];
for i=1:length(ida)
        idx=alldata(ida(i),:)>0.6;
    expression=alldata(ida(i),idx & idx2);
    goodcell=cellID(idx & idx2);
    temp=[];
    for j=1:length(goodcell)
         index=find(strcmp(all1.cell,goodcell(j)));
         temp(j)=celldis(index);
    end
   x1=[sum(temp>0.68), sum(temp>0.55 & temp<=0.68),sum(temp>0.31 & temp<=0.55),sum(temp>0 & temp<=0.31)];
   for j=1:4
       Newdepth(i,j)=x1(j)/total(j);
   end
end
imaging=Newdepth;

%%
[~,~,all]=xlsread('\\rembrandt\Analysis\InSituChipSeq\K4_mm10_P0\MERFISH_comp.xlsx',2);
Name=all(1:end,1);
imaging=cell2mat(all(1:end,2:5));
seq=cell2mat(all(1:end,10:13));
seq2=cell2mat(all(1:end,16:19));
pvalue=zeros(1,length(imaging));
for i=1:length(imaging)
    R=corrcoef(imaging(i,:),seq(i,:));
    pvalue(i)=R(1,2);
end


%%

newcount2=imaging;
for i=1:length(imaging)
    newcount2(i,:)=(imaging(i,:)-mean(imaging(i,:)))./std(imaging(i,:));
end
[~,maxid]=max(newcount2,[],2);
index=1;
name=[];
label=[];
sortedregioncount1=newcount2;
for i=1:4
    idx2=maxid==i;
    [temp,idx3]=sortrows(newcount2(idx2,:),[i],'descend');
    sortedregioncount1(index:index+sum(idx2)-1,:)=temp;
    index=index+sum(idx2);
    tempid=find(idx2);
    name=[name;tempid(idx3)];
    label=[label;Name(tempid(idx3))];
end

figure;
heatmap({'II/III','IV','V','VI'},label,sortedregioncount1,'Colormap',redbluecmap);
grid off;

newcount2=seq;
for i=1:length(seq)
    newcount2(i,:)=(seq(i,:)-mean(seq(i,:)))./std(seq(i,:));
end
figure;
heatmap({'II/III','IV','V','VI'},label,newcount2(name,:),'Colormap',customcolormap_preset('purple-white-green'));
grid off;

pvalue3=pvalue(name);
% newcount3=seq2;
% for i=1:length(seq2)
%     newcount3(i,:)=(seq2(i,:)-mean(seq2(i,:)))./std(seq2(i,:));
% end
% figure;
% heatmap(label,{'II/III','IV','V','VI'},newcount3(name,:)','Colormap',customcolormap_preset('purple-white-green'));
% grid off;
figure;
for i=1:length(pvalue3)
    if pvalue3(i)<0
        sz=5;
    elseif pvalue3(i)<0.3
        sz=15;
    elseif pvalue3(i)<0.8
        sz=30;
    else
        sz=45;
    end
    scatter(1,i,sz,'b','filled');
    hold on;
end
axis off
%%
figure;
heatmap(label,{'II/III','IV','V','VI'},sortedregioncount1','Colormap',redbluecmap);
grid off;
figure;
heatmap(label,{'II/III','IV','V','VI'},newcount2(name,:)','Colormap',customcolormap_preset('purple-white-green'));
grid off;
%%
figure;

