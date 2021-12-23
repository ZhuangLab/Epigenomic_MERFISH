barcodes = GenerateExtendedHammingWords(18);
barcodes = barcodes(sum(barcodes,2)==4,:); 
[~,~,all]=xlsread('\\rembrandt\Analysis\InSituChipSeq\K27ac_mm10_enhancer_E15\Codebook.xlsx',1);
chr=[];
for i=1:length(all)
    temp=strfind(all{i,1},'.');
    chr(i)=str2num(all{i,1}(4:temp(1)));
end
%%
finalbarcodes=[];
finalbad=100;
for i=1:1000000
        barcodes = barcodes(randperm(size(barcodes,1)),:);
        accept=1;
        bad=0;
        for j=1:22
           idx=chr==j;
           tempcode=sum(barcodes(idx,:));
           if  sum(tempcode>5)>0
               accept=0;
               break;
           end
           if sum(tempcode>4)>0
               bad=bad+1;               
           end
        end
        
        if accept && bad<finalbad
            finalbarcodes=barcodes;
            finalbad=bad;
        end
end
%
for j=1:22
           idx=chr==j;
           tempcode=sum(finalbarcodes(idx,:));
           disp( ['Chr:' num2str(j)]);
           disp(tempcode);
           
end


disp(finalbad)