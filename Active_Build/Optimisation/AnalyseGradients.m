function [sortDat]=AnalyseGradients(dirDat,directorySave)
    
    dirDat=ExtractNames(dirDat);
    
    sortDat=GroupData(dirDat);
    for ii=1:numel(sortDat)
        h=GenerateGroupedPlots(sortDat(ii));
        for jj=1:numel(h)
            
            hgsave(h(jj),[directorySave,filesep,h(jj).Name,'.fig'])
        end
    end
    
end

function [h]=GenerateGroupedPlots(sortDat)
    
    diffStep=[sortDat.dirDat(:).diffStep];
    
    diffStep=-(log10(diffStep));
    [diffStepOrd,ordInd]=sort(diffStep);
    
    markList='x+o*sdv^<>ph';
    colorNum=7;
    
    for ii=1:numel(sortDat.dirDat)
        
        gradStart(ordInd(ii),1:numel(sortDat.dirDat(ii).dirChange.grads(1,:)))...
            =sortDat.dirDat(ii).dirChange.grads(1,:);

        gradEnd(ordInd(ii),1:numel(sortDat.dirDat(ii).dirChange.grads(end,:)))...
            =sortDat.dirDat(ii).dirChange.grads(end,:);
        
    end
    h(1)=figure('Name',['GradientData_',sortDat.pattern],'Position',[100 15 800 800]);
    for ii=1:4
        axh(ii)=subplot(2,2,ii);
        hold on
    end
    
    axes(axh(1))
    xlabel('-log10(diffStep)')
    ylabel('Starting Gradient')
    for ii=1:size(gradStart,2)
        l(ii)=plot(diffStepOrd,gradStart(:,ii),['-',markList(mod(ceil(ii/colorNum)-1,13)+1)]);
        l(ii).DisplayName=['var',int2str(ii)];
    end
    legend(l,'Location','EastOutside');
    
    axes(axh(2))
    xlabel('-log10(diffStep)')
    ylabel('Last Gradient')
    for ii=1:size(gradStart,2)
        plot(diffStepOrd,gradEnd(:,ii),['-',markList(mod(ceil(ii/colorNum)-1,13)+1)]);
    end
    
    axes(axh(3))
    xlabel('-log10(diffStep)')
    ylabel('Change in Starting gradient')
    [uniqNum,gradsDiff,meanGradsDiff]=BuildGradDiff(diffStepOrd,gradStart);
    for ii=1:size(gradStart,2)
        l=plot(uniqNum,meanGradsDiff(:,ii)...
            ,['-',markList(mod(ceil(ii/colorNum)-1,13)+1)]); %./abs([1,diffStepOrd(1:end-1)-diffStepOrd(2:end)])'
        l.DisplayName=['var',int2str(ii)];
        l2=plot(diffStepOrd,gradsDiff(:,ii)...
            ,[markList(mod(ceil(ii/colorNum)-1,13)+1)]);
        l2.Color=l.Color;
    end
    
    axes(axh(4))
    xlabel('-log10(diffStep)')
    ylabel('Change in Final gradient')
    [uniqNum,gradsDiff,meanGradsDiff]=BuildGradDiff(diffStepOrd,gradEnd);
    for ii=1:size(gradEnd,2)
        l=plot(uniqNum,meanGradsDiff(:,ii)...
            ,['-',markList(mod(ceil(ii/colorNum)-1,13)+1)]); %./abs([1,diffStepOrd(1:end-1)-diffStepOrd(2:end)])'
        l.DisplayName=['var',int2str(ii)];
        l2=plot(diffStepOrd,gradsDiff(:,ii)...
            ,[markList(mod(ceil(ii/colorNum)-1,13)+1)]);
        l2.Color=l.Color;
    end
    
    h(2)=figure('Name',['GradientEvolution_',sortDat.pattern],'Position',[100 15 800 800]);
    for ii=[1 3 4]
        axh(ii)=subplot(2,2,ii);
        hold on
    end
    axh(2)=subplot(4,2,2);
    hold on
    axh(5)=subplot(4,2,4);
    hold on
    
    axes(axh(1))
    xlabel('Optimisation iteration')
    ylabel('step length')
    for ii=1:numel(diffStepOrd)
        iterNum=numel(sortDat.dirDat(ii).dirChange.Pos);
        l(ii)=plot(1:iterNum,sortDat.dirDat(ii).dirChange.Pos,...
            ['-',markList(mod(ceil(ii/colorNum)-1,13)+1)]);
        l(ii).DisplayName=num2str(-diffStep(ii));
    end
    axh(1).YScale='log';
    legend(l,'Location','NorthEast')
    
    axes(axh(2))
    ylabel('Objective')
    for ii=1:numel(diffStepOrd)
        iterNum=numel(sortDat.dirDat(ii).supportOptim.hist);
        l(ii)=plot(1:iterNum,[sortDat.dirDat(ii).supportOptim.hist(:).obj_curr],...
            ['-',markList(mod(ceil(ii/colorNum)-1,13)+1)]);
        l(ii).DisplayName=num2str(-diffStep(ii));
    end 
    axes(axh(5))
    xlabel('Optimisation iteration')
    ylabel('Gradient Norm')
    for ii=1:numel(diffStepOrd)
        iterNum=numel(sortDat.dirDat(ii).dirChange.gradNorm);
        l(ii)=plot(1:iterNum,[sortDat.dirDat(ii).dirChange.gradNorm],...
            ['-',markList(mod(ceil(ii/colorNum)-1,13)+1)]);
        l(ii).DisplayName=num2str(-diffStep(ii));
    end 
    axh(5).YScale='log';
    axes(axh(3))
    xlabel('Optimisation iteration')
    ylabel('Direction Change - Gradient')
    for ii=1:numel(diffStepOrd)
        iterNum=numel(sortDat.dirDat(ii).dirChange.Grad);
        l(ii)=plot(1:iterNum,sortDat.dirDat(ii).dirChange.Grad,...
            ['-',markList(mod(ceil(ii/colorNum)-1,13)+1)]);
        l(ii).DisplayName=num2str(-diffStep(ii));
    end
    box=axis;
    axis([box(1:2),-1 ,1])
    l3(1)=plot(box(1:2),[1 1],'k-');
    l3(1).DisplayName='Parralel';
    l3(2)=plot(box(1:2),[0 0],'k--');
    l3(2).DisplayName='Normal';
    l3(3)=plot(box(1:2),[-1 -1],'k-.');
    l3(3).DisplayName='Reverse';
    legend(l3,'Location','Southwest')
    
    axes(axh(4))
    xlabel('Optimisation iteration')
    ylabel('Direction Change - Step')
    for ii=1:numel(diffStepOrd)
        iterNum=numel(sortDat.dirDat(ii).dirChange.Step);
        l(ii)=plot(1:iterNum,sortDat.dirDat(ii).dirChange.Step,...
            ['-',markList(mod(ceil(ii/colorNum)-1,13)+1)]);
        l(ii).DisplayName=num2str(-diffStep(ii));
    end
    box=axis;
    axis([box(1:2),-1 ,1])
    l3(1)=plot(box(1:2),[1 1],'k-');
    l3(1).DisplayName='Parralel';
    l3(2)=plot(box(1:2),[0 0],'k--');
    l3(2).DisplayName='Normal';
    l3(3)=plot(box(1:2),[-1 -1],'k-.');
    l3(3).DisplayName='Reverse';
    legend(l3,'Location','Southwest')
end

function [uniqNum,gradsDiff,meanGradsDiff]=BuildGradDiff(diffStep,grads)
    
    diffStepComp=[diffStep;round(diffStep)];
    uniqNum=RemoveIdenticalEntries(diffStepComp(2,:));
    
    for ii=1:numel(uniqNum)
        diffStepComp(3,uniqNum(ii)==diffStepComp(2,:))=ii;
    end
    for ii=1:numel(uniqNum)
        meanGrads(ii,1:size(grads,2))=mean(grads(find(ii==diffStepComp(3,:)),:),1);
    end
    
    meanGradsExp=[zeros([1,size(meanGrads,2)]);meanGrads];
    
%     gradsDiff=log10(abs((grads-meanGradsExp(diffStepComp(3,:),:))./grads));
%     meanGradsDiff=log10(abs((meanGradsExp(1:end-1,:)-meanGradsExp(2:end,:))./meanGradsExp(1:end-1,:)));
    gradsDiff=log10(abs((grads-meanGradsExp(diffStepComp(3,:),:))));
    meanGradsDiff=log10(abs((meanGradsExp(1:end-1,:)-meanGradsExp(2:end,:))));
    
end

function [sortDat]=GroupData(dirDat)
    
    patternList={dirDat(1).pattern};
    patternGroup={1};
    
    for ii=2:numel(dirDat)
        
        flag=false;
        jj=0;
        jjMax=numel(patternList);
        while ~flag && jj<jjMax
            jj=jj+1;
            flag=strcmp(dirDat(ii).pattern,patternList{jj});
            
        end
        
        if flag
           patternGroup{jj}(end+1)=ii;
        else
            patternList{end+1}=dirDat(ii).pattern;
            patternGroup{end+1}=ii;
        end
    end
    
    sortDat=repmat(struct('pattern','','dirDat',struct([])),[1,numel(patternList)]);
    for ii=1:numel(patternList)
        sortDat(ii).pattern=patternList{ii};
        sortDat(ii).dirDat=[dirDat(patternGroup{ii})];
    end
end

function dirDat=ExtractNames(dirDat)
    
    
    
    charNames=char(dirDat(:).case)';
    inds=(regexp(charNames(:)','_'));
    [i,j]=ind2sub(size(charNames),inds);
    keep=true(size(j));
    keep(2:end)=j(1:end-1)~=j(2:end);
    j=j(keep);
    i=i(keep);
    charNames=charNames';
    for ii=1:length(dirDat)
        
        dirDat(ii).pattern=charNames(j(ii),1:i(ii)-1);
        
    end
    
end