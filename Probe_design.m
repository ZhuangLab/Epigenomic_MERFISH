

[~,~,all2]=xlsread('\\rembrandt\Analysis\InSituChipSeq\K27ac_mm10_enhancer_E15\Candidate.xlsx',1);
temp=all2(2:189,1);
TargetStart=[];
TargetEnd=[];
TargetChrom={};
for i=1:length(temp)
    st=temp{i};
    tempst=strsplit(st,'.');
    if length(tempst)==4
        TargetStart(i)=str2num(tempst{3});
    TargetEnd(i)=str2num(tempst{4});
    TargetChrom{i}=tempst{2};
    else
    TargetStart(i)=str2num(tempst{2});
    TargetEnd(i)=str2num(tempst{3});
    TargetChrom{i}=tempst{1};
    
    end
end

allStart=[];
allEnd=[];
allChrom={};

[~,~,all1]=xlsread('\\rembrandt\Analysis\InSituChipSeq\K27ac_mm10_enhancer_E15\Candidate.xlsx',2);
temp=all1(1:end,1);
for i=1:76789
    st=temp{i};
    tempst=strsplit(st,'.');
    allStart(i)=str2num(tempst{2});
    allEnd(i)=str2num(tempst{3});
    allChrom{i}=tempst{1};
end

[~,~,all1]=xlsread('\\rembrandt\Analysis\InSituChipSeq\K27ac_mm10_enhancer_E15\Candidate.xlsx',3);
previous=length(allEnd);
temp=all1(2:end,1);
for i=1:length(temp)
    st=temp{i};
    tempst=strsplit(st,'.');
    allStart(i+previous)=str2num(tempst{3});
    allEnd(i+previous)=str2num(tempst{4});
    allChrom{i+previous}=tempst{2};
end


[~,~,all3]=xlsread('\\rembrandt\Analysis\InSituChipSeq\P300_lib3\Candidates.xlsx',2);
iggStart=cell2mat(all3(1:end,5));
iggEnd=cell2mat(all3(1:end,6));
iggChrom=all3(1:end,4);


genomeseq=fastaread('\\rembrandt\Analysis\InSituChipSeq\P300lib2\mm10.fa');
savePath='\\rembrandt\Analysis\InSituChipSeq\K27ac_mm10_enhancer_E15\';
%%
Target={};
for i=1:length(TargetStart)
    temp=genomeseq(find(strcmp({genomeseq.Header},TargetChrom{i}))).Sequence;
    Target{i}=temp(TargetStart(i):TargetEnd(i));
end

allseq={};
for i=1:length(allStart)
    temp=genomeseq(find(strcmp({genomeseq.Header},allChrom{i}))).Sequence;
    allseq{i}=temp(allStart(i):allEnd(i));
end
%%
iggseq={};
for i=1:length(iggStart)
    temp=genomeseq(find(strcmp({genomeseq.Header},iggChrom{i}))).Sequence;
    iggseq{i}=temp(iggStart(i):iggEnd(i));
end

%%  design probes for each target
warning('off');
GClow=33;
GChigh=73;
TMlow=61;
TMhigh=81;
similarthreshold=15;
targetregion=cell(length(TargetStart),1);
regionnum=zeros(1,length(TargetStart));
targetregionstart=cell(length(TargetStart),1);
targetregionend=cell(length(TargetStart),1);
parfor i=1:length(TargetStart)
    temp=Target{i};
    index=1;
    k=1;
    temps={};
    temp1=[];
    temp2=[];
    strand=1;
    while k<length(temp)-30    
            candidate=temp(k:k+29);   
            GC=(length(strfind(candidate,'G'))+length(strfind(candidate,'C')))/length(candidate)*100;
            TM=64.9 +41*(length(strfind(candidate,'G'))+length(strfind(candidate,'C'))-16.4)/length(candidate);
            if GC>=GClow && GC<=GChigh && TM>=TMlow && TM<=TMhigh && ~contains(candidate,{'AAAA','GGGG','CCCC','TTTT','AGAGAG','ACACAC','TCTCTC','TGTGTG','ATATAT'})
            check=1;   
            j=1;        
                    while j<35000
                        m=allseq{j};
                        if abs(allStart(j)-TargetStart(i))>1000 | strcmp(allChrom{j},TargetChrom{i})==0
                            [score,align]=swalign(candidate,m);       
                            if max(diff(find([1,diff(align(2,:)),1])))>similarthreshold  % longest consecutive match
                                check=0;
                                break;
                            end
                            [score,align]=swalign(candidate,seqrcomplement(m));       
                            if max(diff(find([1,diff(align(2,:)),1])))>similarthreshold  % longest consecutive match
                                check=0;
                                break;
                            end
                        end
                        j=j+1;
                    end
                 j=1;
                    while j<7500
                        m=iggseq{j};
                        
                        [score,align]=swalign(candidate,m);       
                        if max(diff(find([1,diff(align(2,:)),1])))>similarthreshold  % longest consecutive match
                            check=0;
                            break;
                        end
                        [score,align]=swalign(candidate,seqrcomplement(m));       
                        if max(diff(find([1,diff(align(2,:)),1])))>similarthreshold  % longest consecutive match
                            check=0;
                            break;
                        end
                        
                        j=j+1;
                    end 
%                     disp(j);
                    if check
                        temps{index}=candidate;
                        temp1(index)=k;
                        temp2(index)=k+29;
                        index=index+1;
                        k=k+27;
%                         disp(candidate);
%                         disp(GC);
%                         disp(TM);
                    
                    end
               
            end
            k=k+1;
    end
                targetregion{i}=temps;
            targetregionstart{i}=temp1;
            targetregionend{i}=temp2;
    regionnum(i)=index-1;
    disp(['Find ' num2str(index-1) ' probes for region ' num2str(i) ' from length ' num2str(length(temp)) ' nt']);
end
save([savePath '\targetregion.mat'],'targetregion','regionnum');
%% quantify probe number
figHandle = figure('Name', 'Number of probes per loci', 'Color', 'w');

[~, sind] = sort(regionnum, 'Descend');
plot(regionnum(sind), '.'); hold on;
xlabel('loci ID');
ylabel('Number of regions');
SaveFigure(figHandle, 'overwrite', true, 'formats', {'fig', 'png'},'savePath', savePath);

%% quantify probe number per 100 bp
figHandle = figure('Name', 'Number of probes per 100bp', 'Color', 'w');
temp=regionnum./(TargetEnd-TargetStart)*100;
[~, sind] = sort(temp, 'Descend');
plot(temp(sind), '.'); hold on;

xlabel('loci ID');
ylabel('Number of probes per 100 bp');
title(['Avg number of probes per 100 bp:', num2str(mean(regionnum'./(TargetEnd-TargetStart))*100)]);
SaveFigure(figHandle, 'overwrite', true, 'formats', {'fig', 'png'},'savePath', savePath);


%% Build oligos


libraryName='K27ac_mm10_E15';
[~,~,all5]=xlsread('\\rembrandt\Analysis\InSituChipSeq\H3K27ac_660loci_RPE\readout.xlsx',1);
readouts = all5(2:end,3);
readoutname=all5(2:end,2);

[~,~,all4]=xlsread('\\rembrandt\Analysis\InSituChipSeq\K27ac_mm10_enhancer_E15\Codebook',1);
temp=all4(1:147,1);
newbarcodes=cell2mat(all4(1:147,2:end));
for i=1:length(temp)
    st=temp{i};
    tempst=strsplit(st,'.');
    localGenestart(i)=str2num(tempst{2});
    localGeneend(i)=str2num(tempst{3});
    localGenechr{i}=tempst{1};
end

oligosPath = [savePath 'oligos.fasta'];
% if ~exist(oligosPath)
    oligos = [];

for i=1:length(localGenechr)
        % Display progress
        localGeneName=[localGenechr{i} ':' num2str(localGenestart(i)) '-' num2str(localGeneend(i))];
        PageBreak();
        display(['Designing probes for ' libraryName ': ' localGeneName]);
        ids2=find(TargetEnd==localGeneend(i) & TargetStart==localGenestart(i) & strcmp(TargetChrom,localGenechr{i}));
        ids2=ids2(1);
        % Determine the bits to include for each word
        possibleReadouts = readouts(newbarcodes(i,:)==1); 
        readName=readoutname(newbarcodes(i,:)==1);
            seqs = {};
            headers = {};
            missing=4;
            % Build all possible oligos
            for p=1:regionnum(ids2)
                % Create random orientation and selection of readouts
                id=boolean([1,1,1,1]);
                id(missing)=0;
                localReadouts = possibleReadouts(id);
                seq=seqrcomplement(targetregion{ids2}{p});
                if mod(p,2)==0
                    seq=targetregion{ids2}{p};
                end
                if rand(1) > 0.5
                    % Create header 
                    headers{p} = [libraryName ' ' ...
                        'A ' ... % Add A pad
                        readName{1} ' ' ...
                        readName{2} ' '...
                        localGeneName ' ' ...
                        'A ' ... % Add A pad
                        readName{3}];

                    % Create sequence
                    seqs{p} = ['A ' seqrcomplement(localReadouts{1}) ' '...
                           'A ' seqrcomplement(localReadouts{2}) ' A '...
                            seq ' A ' ...
                            seqrcomplement(localReadouts{3})];
                else
                   headers{p} = [libraryName ' ' ...
                        'A ' ... % Add A pad
                        readName{3} ' ' ...
                        localGeneName ' ' ...
                        'A ' ... % Add A pad
                        readName{1} ' ' ...
                        readName{2}];

                    % Create sequence
                    seqs{p} = ['A ' seqrcomplement(localReadouts{3}) ' A '... 
                            seq ' A ' ...
                            seqrcomplement(localReadouts{1}) ' A ' ...
                            seqrcomplement(localReadouts{2})];
                end
                missing=missing-1;
                if missing==0
                    missing=4;
                end
            end
            display(['... constructed ' num2str(length(seqs)) ' possible probes']);


            % Save new oligos in oligos struct
            for s=1:length(seqs)
                oligos(end+1).Header = headers{s};
                oligos(end).Sequence = seqs{s};
            end
 end
    PageBreak();
    display(['Writing: ' oligosPath]);
    writeTimer = tic;
    fastawrite(oligosPath, oligos);
    display(['... completed in ' num2str(toc(writeTimer))]);


%% Load fasta file and select primers for library

primersPath='\\rembrandt\Analysis\InSituChipSeq\H3K27me3_lib\H3K27me3_primers.fasta';
primers = fastaread(primersPath);
usedPrimers = primers(1:2);

% Add primers to encoding probes
PageBreak();
display('Adding primers');
finalPrimersPath = [savePath libraryName '_primers.fasta'];
% if ~exist(finalPrimersPath)
    % Record the used primers
    fastawrite(finalPrimersPath, usedPrimers);
    display(['Wrote: ' finalPrimersPath]);

    % Build the final oligos
    finalOligos = [];
    for i=1:length(oligos)
        stringParts = strsplit(oligos(i).Header, ' ');
        name1 = strsplit(usedPrimers(1).Header, ' ');
        name1 = name1{1};
        name2 = strsplit(usedPrimers(2).Header, ' ');
        name2 = name2{1};
        finalOligos(i).Header = [stringParts{1} ' ' ...
        name1 ' '];
        for j=2:length(stringParts)
            finalOligos(i).Header = [finalOligos(i).Header ...
                stringParts{j} ' '];
        end
        finalOligos(i).Header = [finalOligos(i).Header name2];

        finalOligos(i).Sequence = [usedPrimers(1).Sequence ' ' ...
        oligos(i).Sequence ' ' ...
        seqrcomplement(usedPrimers(2).Sequence)];
    end
    
    
%% Final cross checks and writing final file
PageBreak();
rRNAOtTable = OTTable.Load('E:\Users\Won\lib_construction_mouse\Mus_musculus_UCSC_mm10\mm10_wj.kallisto_out\v4_2_analysis\TRDesigner_analysis\trDesigner\OTTable_1');
display('Running final cross checks and building final fasta file');
% Write final fasta
oligosPath = [savePath libraryName '_oligos.fasta'];
% if ~exist(oligosPath)
    % Screen against the original tables
    tic;
    display(['Searching oligos for homology']);
    hasrRNAPenalty = cellfun(@(x) sum(rRNAOtTable.CalculatePenalty(seqrcomplement(x(~isspace(x)))))>0, {finalOligos.Sequence});
    display(['... completed in ' num2str(toc) ' s']);

    indsToKeep = ~hasrRNAPenalty;
    display(['... found ' num2str(sum(~indsToKeep)) ' oligos to remove ']);
    indsToRemove = find(~indsToKeep);
    for r=1:length(indsToRemove)
        display(['...     ' finalOligos(indsToRemove(r)).Header]);
    end

    % Remove bad oligos
    finalOligos = finalOligos(indsToKeep);

    % Write final oligos
    fastawrite(oligosPath, finalOligos);
    display(['Wrote: ' oligosPath]);