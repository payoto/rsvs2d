
function [h]=OptimHistory(isGradient,optimstruct,knownOptim,defaultVal,dirOptim)
    
    if ~isGradient
      
            [h]=OptimHistory_nograd(optimstruct,knownOptim,defaultVal,dirOptim);
    else
            [h]=OptimHistory_grad(optimstruct,knownOptim,defaultVal,dirOptim);
    end
    
end
%{
function [h]=OptimHistory_old(optimstruct,knownOptim,dirOptim)
    
    
    h=figure('Name','Optimisation Result','Position',[20 100 1000 600]);
    
    % Plot 1
    subplot(1,2,1,'ticklabelinterpreter','latex')
    nVar=length(optimstruct(1).population);
    nIter=length(optimstruct);
    iterRes=zeros([nIter,nVar]);
    hold on
    for ii=1:nIter
        nVarLoc=length(optimstruct(ii).population);
        iterRes(ii,1:nVarLoc)=[optimstruct(ii).population(:).objective];
        lSub1(1)=plot(ones(1,nVarLoc)*ii,iterRes(ii,1:nVarLoc),'b.','markersize',5);
    end
    switch dirOptim
        case 'min'
            iterRes(iterRes==0)=1000;
            minRes=min(iterRes,[],2);
        case 'max'
            iterRes(iterRes==0)=-1000;
            minRes=max(iterRes,[],2);
    end
    
    meanRes=mean(iterRes,2);
    stdRes=std(iterRes,0,2);
    lSub1(2)=plot(1:nIter,minRes,'r-');
    lSub1(3)=plot(1:nIter,meanRes,'color',[0.7 0 0]);
    lSub1(4)=plot([0,nIter],[knownOptim knownOptim],'r--');
    
    legend(lSub1,{'Population',['Population ',dirOptim,'imum'],'Population mean',...
        'Theoretical Optimum'},'Location','NorthEast', 'interpreter','latex');
    
    xlabel('Iteration', 'interpreter','latex','fontsize',12)
    ylabel('$J(\mathbf{x})$', 'interpreter','latex','fontsize',12)
    
    switch dirOptim
        case 'min'
            testOrder=max(minRes);
            orderSol=ceil(-log10(abs(testOrder)));
            box(4)=ceil(testOrder*10^orderSol)*10^(-orderSol);
            testOrder=min([minRes;knownOptim]);
            box(3)=floor(testOrder*10^orderSol)*10^(-orderSol);
        case 'max'
            testOrder=min(minRes);
            orderSol=ceil(-log10(abs(testOrder)));
            box(3)=floor(testOrder*10^orderSol)*10^(-orderSol);
            testOrder=max([minRes;knownOptim]);
            box(4)=ceil(testOrder*10^orderSol)*10^(-orderSol);
    end
    box(1:2)=[0,nIter+1];
    
    
    axis(box);
    
    xT=box(2)-(box(2)-box(1))*0.05;
    yT=min(minRes)-(box(4)-box(3))*0.05;
    strT=['$\quad\quad$ $J^*(\mathbf{x})$ = ',sprintf('%10.3e',(min(minRes)))];
    strT={strT,['$J^*(\mathbf{x})-J^*_T$ = ',sprintf('%10.3e',min(minRes)-knownOptim)]};

    strT=regexprep(strT,'\ ','\\space');
    text(xT,yT,strT, 'interpreter','latex','HorizontalAlignment','right');
    
    
    % Plot 2
    axh=subplot(1,2,2,'ticklabelinterpreter','latex');
    
    switch dirOptim
        case 'min'
            errMeasure=-(knownOptim-minRes);
            errMean=-(knownOptim-meanRes);
        case 'max'
            errMeasure=knownOptim-minRes;
            errMean=(knownOptim-meanRes);
    end
    
    lSub2(1)=semilogy(1:nIter,errMeasure);
    hold on
    lSub2(2)=semilogy(1:nIter,errMean);
    
    for ii=1:nIter
        stdline=[errMean(ii)-stdRes(ii),errMean(ii)+stdRes(ii)];
        lSub2(3)=semilogy([ii ii],stdline,'g:+');
    end
    legend(lSub2,{['Population ',dirOptim,'imum'],'Population mean',...
        'Standard Deviation'},'Location','NorthEast', 'interpreter','latex');
    
    xlabel('Iteration', 'interpreter','latex','fontsize',12)
    ylabel('$J^*_T-J(\mathbf{x})$', 'interpreter','latex','fontsize',12)
    set(axh,'ticklabelinterpreter','latex')
    
end

function [h]=OptimHistory_nograd(optimstruct,knownOptim,defaultVal,dirOptim)
    
    
    h=figure('Name','Optimisation Result','Position',[20 100 1000 600]);
    
    % Plot 1
    subplot(1,2,1,'ticklabelinterpreter','latex')
    nVar=length(optimstruct(1).population);
    nIter=length(optimstruct);
    iterRes=zeros([nIter,nVar])+defaultVal;
    hold on
    
    for ii=1:nIter
        nVarLoc=length(optimstruct(ii).population);
        iterRes(ii,1:nVarLoc)=[optimstruct(ii).population(:).objective];
        lSub1(1)=plot(ones(1,nVarLoc)*ii,iterRes(ii,1:nVarLoc),'b.','markersize',5);
    end
    switch dirOptim
        case 'min'
            [minRes,minPos]=min(iterRes,[],2);
        case 'max'
            [minRes,minPos]=max(iterRes,[],2);
    end
    
    meanRes=mean(iterRes,2);
    stdRes=std(iterRes,0,2);
    lSub1(2)=plot(1:nIter,minRes,'r-');
    lSub1(3)=plot(1:nIter,meanRes,'color',[0.7 0 0]);
    lSub1(4)=plot([0,nIter],[knownOptim knownOptim],'r--');
    
    legend(lSub1,{'Population',['Population ',dirOptim,'imum'],'Population mean',...
        'Theoretical Optimum'},'Location','NorthEast', 'interpreter','latex');
    
    xlabel('Iteration', 'interpreter','latex','fontsize',12)
    ylabel('$J(\mathbf{x})$', 'interpreter','latex','fontsize',12)
    
    switch dirOptim
        case 'min'
            testOrder=max(minRes);
            orderSol=ceil(-log10(abs(testOrder)));
            box(4)=ceil(testOrder*10^orderSol)*10^(-orderSol);
            testOrder=min([minRes;knownOptim]);
            box(3)=floor(testOrder*10^orderSol)*10^(-orderSol);
        case 'max'
            testOrder=min(minRes);
            orderSol=ceil(-log10(abs(testOrder)));
            box(3)=floor(testOrder*10^orderSol)*10^(-orderSol);
            testOrder=max([minRes;knownOptim]);
            box(4)=ceil(testOrder*10^orderSol)*10^(-orderSol);
    end
    box(1:2)=[0,nIter+1];
    
    
    axis(box);
    
    xT=box(2)-(box(2)-box(1))*0.05;
    yT=min(minRes)-(box(4)-box(3))*0.05;
    strT=['$\quad\quad$ $J^*(\mathbf{x})$ = ',sprintf('%10.3e',(min(minRes)))];
    strT={strT,['$J^*(\mathbf{x})-J^*_T$ = ',sprintf('%10.3e',min(minRes)-knownOptim)]};

    strT=regexprep(strT,'\ ','\\space');
    text(xT,yT,strT, 'interpreter','latex','HorizontalAlignment','right');
    
    
    % Plot 2
    
    % Plot 2
    axh=subplot(1,2,2,'ticklabelinterpreter','latex');
    
    switch dirOptim
        case 'min'
            errMeasure=-(knownOptim-minRes);
            
        case 'max'
            errMeasure=knownOptim-minRes;
    end

    
    popNominal=zeros([nIter,numel(optimstruct(1).population(1).fill)]);
    objNominal=zeros([nIter,1]);
    cBounds=[min(minRes),max(minRes)];
    [datCol]=ProjectColormap(h.Colormap,meanRes,cBounds);
    
    for ii=1:nIter
        nPop=length(optimstruct(ii).population);
        for jj=1:nPop
            plot3(ones(size(optimstruct(ii).population(jj).fill))*((ii+1)/2),...
                1:numel(optimstruct(ii).population(jj).fill),...
                optimstruct(ii).population(jj).fill,'Color',datCol(ii,:));
            hold on
        end
    end
    
    for ii=1:nIter
        
        popNominal(ii,:)=optimstruct(ii).population(minPos(ii)).fill;
        objNominal(ii)=optimstruct(ii).population(minPos(ii)).objective;
    end
    
    
    iterVec=1:nIter;
    varVec=1:length(popNominal(1,:));
    
    [iterGrid,varGrid]=meshgrid(iterVec,varVec);
    [objGrid,~]=meshgrid(objNominal,varVec);
    
    surf(iterGrid,varGrid,popNominal',objGrid)
    
    colorbar;
    
    xlabel('Iteration', 'interpreter','latex','fontsize',12)
    ylabel('Variable Index', 'interpreter','latex','fontsize',12)
    zlabel('Variable Value', 'interpreter','latex','fontsize',12)
    set(axh,'ticklabelinterpreter','latex')
    
    %{
    axh=subplot(1,2,2,'ticklabelinterpreter','latex');
    
    switch dirOptim
        case 'min'
            errMeasure=-(knownOptim-minRes);
            errMean=-(knownOptim-meanRes);
        case 'max'
            errMeasure=knownOptim-minRes;
            errMean=(knownOptim-meanRes);
    end
    
    lSub2(1)=semilogy(1:nIter,errMeasure);
    hold on
    lSub2(2)=semilogy(1:nIter,errMean);
    
    for ii=1:nIter
        stdline=[errMean(ii)-stdRes(ii),errMean(ii)+stdRes(ii)];
        lSub2(3)=semilogy([ii ii],stdline,'g:+');
    end
    legend(lSub2,{['Population ',dirOptim,'imum'],'Population mean',...
        'Standard Deviation'},'Location','NorthEast', 'interpreter','latex');
    
    xlabel('Iteration', 'interpreter','latex','fontsize',12)
    ylabel('$J^*_T-J(\mathbf{x})$', 'interpreter','latex','fontsize',12)
    set(axh,'ticklabelinterpreter','latex')
    %}
    
    
end

function [h]=OptimHistory_grad(optimstruct,knownOptim,defaultVal,dirOptim)
    
    h=figure('Name','Optimisation Result','Position',[20 100 1000 600]);
    
    % Plot 1
    subplot(1,2,1,'ticklabelinterpreter','latex')
    nVar=length(optimstruct(1).population);
    nIter=length(optimstruct);
    iterRes=zeros([nIter,nVar])+defaultVal;
    hold on
    for ii=1:nIter
        nVarLoc=length(optimstruct(ii).population);
        iterRes(ii,1:nVarLoc)=[optimstruct(ii).population(:).objective];
        %lSub1(1)=plot(ones(1,nVarLoc)*ii,iterRes(ii,1:nVarLoc),'b.','markersize',5);
    end
    for ii=2:2:nIter
        nVarLoc=length(optimstruct(ii).population);
        lSub1(1)=plot(ones(1,nVarLoc)*(ii+2)/2,iterRes(ii,1:nVarLoc),'b.','markersize',5);
    end
    switch dirOptim
        case 'min'
            minRes=iterRes(1:2:end,1);
        case 'max'
            minRes=iterRes(1:2:end,1);
    end
    
    meanRes=mean(iterRes,2);
    stdRes=std(iterRes,0,2);
    lSub1(2)=plot([2:2:nIter]/2,minRes,'r-');
    %lSub1(3)=plot(1:2:nIter,meanRes,'color',[0.7 0 0]);
    lSub1(3)=plot([0,nIter],[knownOptim knownOptim],'r--');
    
    legend(lSub1,{'Population',['Root Member'],...
        'Theoretical Optimum'},'Location','NorthEast', 'interpreter','latex');
    
    xlabel('Iteration', 'interpreter','latex','fontsize',12)
    ylabel('$J(\mathbf{x})$', 'interpreter','latex','fontsize',12)
    
    switch dirOptim
        case 'min'
            testOrder=max(minRes);
            orderSol=ceil(-log10(abs(testOrder)));
            box(4)=ceil(testOrder*10^orderSol)*10^(-orderSol);
            testOrder=min([minRes;knownOptim]);
            box(3)=floor(testOrder*10^orderSol)*10^(-orderSol);
        case 'max'
            testOrder=min(minRes);
            orderSol=ceil(-log10(abs(testOrder)));
            box(3)=floor(testOrder*10^orderSol)*10^(-orderSol);
            testOrder=max([minRes;knownOptim]);
            box(4)=ceil(testOrder*10^orderSol)*10^(-orderSol);
    end
    box(1:2)=[0,ceil(nIter/2)+1];
    
    
    axis(box);
    
    xT=box(2)-(box(2)-box(1))*0.05;
    yT=min(minRes)-(box(4)-box(3))*0.05;
    strT=['$\quad\quad$ $J^*(\mathbf{x})$ = ',sprintf('%10.3e',(min(minRes)))];
    strT={strT,['$J^*(\mathbf{x})-J^*_T$ = ',sprintf('%10.3e',min(minRes)-knownOptim)]};

    strT=regexprep(strT,'\ ','\\space');
    text(xT,yT,strT, 'interpreter','latex','HorizontalAlignment','right');
    
    
    % Plot 2
    axh=subplot(1,2,2,'ticklabelinterpreter','latex');
    
    switch dirOptim
        case 'min'
            errMeasure=-(knownOptim-minRes);
        case 'max'
            errMeasure=knownOptim-minRes;
    end
    
    
%     for ii=1:nIter
%         nPop=length(optimstruct(ii).population);
%         for jj=1:nPop
%             plot3(ones(size(optimstruct(ii).population(jj).fill))*((ii+1)/2),...
%                 1:numel(optimstruct(ii).population(jj).fill),...
%                 optimstruct(ii).population(jj).fill)
%             hold on
%         end
%     end
    
    popNominal=zeros([ceil(nIter/2),numel(optimstruct(1).population(1).fill)]);
    objNominal=zeros([ceil(nIter/2),1]);
    kk=0;
    for ii=1:2:nIter
        kk=kk+1;
        popNominal(kk,:)=optimstruct(ii).population(1).fill;
        objNominal(kk)=optimstruct(ii).population(1).objective;
    end
    
    
    iterVec=((2:2:length(optimstruct))-1)/2;
    varVec=1:length(popNominal(1,:));
    
    [iterGrid,varGrid]=meshgrid(iterVec,varVec);
    [objGrid,~]=meshgrid(objNominal,varVec);
    
    surf(iterGrid,varGrid,popNominal',objGrid)
    
    colorbar;
    
    xlabel('Iteration', 'interpreter','latex','fontsize',12)
    ylabel('Variable Index', 'interpreter','latex','fontsize',12)
    zlabel('Variable Value', 'interpreter','latex','fontsize',12)
    set(axh,'ticklabelinterpreter','latex')
    
end

%}


function [h]=OptimHistory_nograd(optimstruct,knownOptim,defaultVal,dirOptim)
    
    
    h(1)=figure('Name','Optimisation Result','Position',[20 100 1000 600]);
    
    % Plot 1
    subplot(1,2,1,'ticklabelinterpreter','latex')
    nVar=length(optimstruct(1).population);
    nIter=length(optimstruct);
    iterRes=zeros([nIter,nVar])+defaultVal;
    hold on
    cVec=[]; 
    for ii=1:length(optimstruct),
        for jj=1:length(optimstruct(ii).population), 
            cVec=[cVec,optimstruct(ii).population(jj).additional.c]; 
        end  
    end 
    
    meanC=mean(cVec);
    for ii=1:nIter
        nVarLoc=length(optimstruct(ii).population);
        iterRes(ii,1:nVarLoc)=[optimstruct(ii).population(:).objective]/meanC;
        
        lSub1(1)=plot(ones(1,nVarLoc)*ii,iterRes(ii,1:nVarLoc),'b.','markersize',5);
    end
    switch dirOptim
        case 'min'
            [minRes,minPos]=min(iterRes,[],2);
        case 'max'
            [minRes,minPos]=max(iterRes,[],2);
    end
    
    meanRes=mean(iterRes,2);
    stdRes=std(iterRes,0,2);
    lSub1(2)=plot(1:nIter,minRes,'r-');
    lSub1(3)=plot(1:nIter,meanRes,'color',[0.7 0 0]);
    lSub1(4)=plot([0,nIter],[knownOptim knownOptim],'r--');
    
    legend(lSub1,{'Population',['Population ',dirOptim,'imum'],'Population mean',...
        'Theoretical Optimum'},'Location','NorthEast', 'interpreter','latex');
    
    xlabel('Iteration', 'interpreter','latex','fontsize',12)
    ylabel('$J(\mathbf{x})$', 'interpreter','latex','fontsize',12)
    
    switch dirOptim
        case 'min'
            testOrder=max(minRes);
            orderSol=ceil(-log10(abs(testOrder)));
            box(4)=ceil(testOrder*10^orderSol)*10^(-orderSol);
            testOrder=min([minRes;knownOptim]);
            box(3)=floor(testOrder*10^orderSol)*10^(-orderSol);
        case 'max'
            testOrder=min(minRes);
            orderSol=ceil(-log10(abs(testOrder)));
            box(3)=floor(testOrder*10^orderSol)*10^(-orderSol);
            testOrder=max([minRes;knownOptim]);
            box(4)=ceil(testOrder*10^orderSol)*10^(-orderSol);
    end
    box(1:2)=[0,nIter+1];
    
    
    axis(box);
    
    xT=box(2)-(box(2)-box(1))*0.05;
    yT=min(minRes)-(box(4)-box(3))*0.05;
    strT=['$\quad\quad$ $J^*(\mathbf{x})$ = ',sprintf('%10.3e',(min(minRes)))];
    strT={strT,['$J^*(\mathbf{x})-J^*_T$ = ',sprintf('%10.3e',min(minRes)-knownOptim)]};

    strT=regexprep(strT,'\ ','\\space');
    text(xT,yT,strT, 'interpreter','latex','HorizontalAlignment','right');
    
    
    % Plot 2
    
    % Plot 2
    axh=subplot(1,2,2,'ticklabelinterpreter','latex');
    
    switch dirOptim
        case 'min'
            errMeasure=-(knownOptim-minRes);
            
        case 'max'
            errMeasure=knownOptim-minRes;
    end

    
    popNominal=zeros([nIter,numel(optimstruct(1).population(1).fill)]);
    objNominal=zeros([nIter,1]);
    cBounds=[min(minRes),max(minRes)];
    [datCol]=ProjectColormap(h.Colormap,meanRes,cBounds);
    maxFill=popNominal;
    minFill=popNominal;
    meanFill=popNominal;
    stdFill=popNominal;
        %             plot3(ones(size(optimstruct(ii).population(jj).fill))*((ii+1)/2),...
%                 1:numel(optimstruct(ii).population(jj).fill),...
%                 optimstruct(ii).population(:).fill,'Color',datCol(ii,:));
     for ii=1:nIter
        
        popNominal(ii,:)=optimstruct(ii).population(minPos(ii)).fill;
        objNominal(ii)=optimstruct(ii).population(minPos(ii)).objective;
        fillPop=vertcat(optimstruct(ii).population(:).fill);
        maxFill(ii,:)=max(fillPop);
        minFill(ii,:)=min(fillPop);
        meanFill(ii,:)=mean(fillPop);
        stdFill(ii,:)=std(fillPop);
        
        plot3(ones(size(optimstruct(ii).population(1).fill))*((ii+1)),...
            1:numel(optimstruct(ii).population(1).fill),...
            maxFill(ii,:),'Color',datCol(ii,:));
        hold on
        plot3(ones(size(optimstruct(ii).population(1).fill))*((ii+1)),...
            1:numel(optimstruct(ii).population(1).fill),...
            minFill(ii,:),'Color',datCol(ii,:));
    end
    
    
    iterVec=1:nIter;
    varVec=1:length(popNominal(1,:));
    
    [iterGrid,varGrid]=meshgrid(iterVec,varVec);
    [objGrid,~]=meshgrid(objNominal,varVec);
    [meanGrid,~]=meshgrid(meanRes,varVec);
    
    surf(iterGrid,varGrid,popNominal',objGrid)
    colorbar;
    xlabel('Iteration', 'interpreter','latex','fontsize',12)
    ylabel('Variable Index', 'interpreter','latex','fontsize',12)
    zlabel('Variable Value', 'interpreter','latex','fontsize',12)
    set(axh,'ticklabelinterpreter','latex')
%     axh=subplot(2,2,4,'ticklabelinterpreter','latex');
%     hold on
%     
%     mesh(iterGrid,varGrid,maxFill')
%     hold on
%     surf(iterGrid,varGrid,minFill')

    % Figure 2
    h(2)=figure('Name','Design Variable Evolution','Position',[20 100 1000 600]);
    
    axh=subplot(2,2,1,'ticklabelinterpreter','latex');
    s(1)=surf(iterGrid,varGrid,meanFill');
    caxis([0 1])
    
    axc(1)=colorbar;
    view(0 ,90)
    xlabel('Iteration', 'interpreter','latex','fontsize',12)
    ylabel('Variable Index', 'interpreter','latex','fontsize',12)
    
    set(axh,'ticklabelinterpreter','latex')
    axh=subplot(2,2,3,'ticklabelinterpreter','latex');
    s(2)=surf(iterGrid,varGrid,stdFill');
    caxis([0 1])
    axc(2)=colorbar;
    view(0 ,90)
    xlabel('Iteration', 'interpreter','latex','fontsize',12)
    ylabel('Variable Index', 'interpreter','latex','fontsize',12)
    
    set(axh,'ticklabelinterpreter','latex')
    axh=subplot(2,2,2,'ticklabelinterpreter','latex');
    s(3)=surf(iterGrid,varGrid,maxFill');
    caxis([0 1])
    axc(3)=colorbar;
    view(0 ,90)
    xlabel('Iteration', 'interpreter','latex','fontsize',12)
    ylabel('Variable Index', 'interpreter','latex','fontsize',12)
    
    axh=subplot(2,2,4,'ticklabelinterpreter','latex');
    s(4)=surf(iterGrid,varGrid,minFill');
    caxis([0 1])
    axc(4)=colorbar;
    view(0 ,90)
    xlabel('Iteration', 'interpreter','latex','fontsize',12)
    ylabel('Variable Index', 'interpreter','latex','fontsize',12)
    
    cLabel={'Mean','Standard Deviation','Maximum','Minimum'};
    for ii=1:length(s)
        s(ii).EdgeColor='none';
        axc(ii).TickLabelInterpreter='latex';
        axc(ii).Label.String=cLabel{ii};
        axc(ii).Label.FontSize=14;
    end
    
    
    
    %{
    axh=subplot(1,2,2,'ticklabelinterpreter','latex');
    
    switch dirOptim
        case 'min'
            errMeasure=-(knownOptim-minRes);
            errMean=-(knownOptim-meanRes);
        case 'max'
            errMeasure=knownOptim-minRes;
            errMean=(knownOptim-meanRes);
    end
    
    lSub2(1)=semilogy(1:nIter,errMeasure);
    hold on
    lSub2(2)=semilogy(1:nIter,errMean);
    
    for ii=1:nIter
        stdline=[errMean(ii)-stdRes(ii),errMean(ii)+stdRes(ii)];
        lSub2(3)=semilogy([ii ii],stdline,'g:+');
    end
    legend(lSub2,{['Population ',dirOptim,'imum'],'Population mean',...
        'Standard Deviation'},'Location','NorthEast', 'interpreter','latex');
    
    xlabel('Iteration', 'interpreter','latex','fontsize',12)
    ylabel('$J^*_T-J(\mathbf{x})$', 'interpreter','latex','fontsize',12)
    set(axh,'ticklabelinterpreter','latex')
    %}
    
    
end

function [h]=OptimHistory_grad(optimstruct,knownOptim,defaultVal,dirOptim)
    
    h=figure('Name','Optimisation Result','Position',[20 100 1000 600]);
    
    % Plot 1
    subplot(1,2,1,'ticklabelinterpreter','latex')
    [iterRes,nIter,nVar]=BuildIterRes(optimstruct,defaultVal);
    hold on
    for ii=1:nIter
        nVarLoc=length(optimstruct(ii).population);
        iterRes(ii,1:nVarLoc)=[optimstruct(ii).population(:).objective];
        %lSub1(1)=plot(ones(1,nVarLoc)*ii,iterRes(ii,1:nVarLoc),'b.','markersize',5);
    end
    for ii=2:2:nIter
        nVarLoc=length(optimstruct(ii).population);
        lSub1(1)=plot(ones(1,nVarLoc)*(ii+2)/2,iterRes(ii,1:nVarLoc),'b.','markersize',5);
    end
    switch dirOptim
        case 'min'
            minRes=iterRes(1:2:end,1);
        case 'max'
            minRes=iterRes(1:2:end,1);
    end
    
    meanRes=mean(iterRes,2);
    stdRes=std(iterRes,0,2);
    lSub1(2)=plot(1:numel(minRes),minRes,'r-');
    %lSub1(3)=plot(1:2:nIter,meanRes,'color',[0.7 0 0]);
    lSub1(3)=plot([0,nIter],[knownOptim knownOptim],'r--');
    
    legend(lSub1,{'Population',['Root Member'],...
        'Theoretical Optimum'},'Location','NorthEast', 'interpreter','latex');
    
    xlabel('Iteration', 'interpreter','latex','fontsize',12)
    ylabel('$J(\mathbf{x})$', 'interpreter','latex','fontsize',12)
    
    switch dirOptim
        case 'min'
            testOrder=max(minRes);
            orderSol=ceil(-log10(abs(testOrder)));
            box(4)=ceil(testOrder*10^orderSol)*10^(-orderSol);
            testOrder=min([minRes;knownOptim]);
            box(3)=floor(testOrder*10^orderSol)*10^(-orderSol);
        case 'max'
            testOrder=min(minRes);
            orderSol=ceil(-log10(abs(testOrder)));
            box(3)=floor(testOrder*10^orderSol)*10^(-orderSol);
            testOrder=max([minRes;knownOptim]);
            box(4)=ceil(testOrder*10^orderSol)*10^(-orderSol);
    end
    box(1:2)=[0,ceil(nIter/2)+1];
    
    
    axis(box);
    
    xT=box(2)-(box(2)-box(1))*0.05;
    yT=min(minRes)-(box(4)-box(3))*0.05;
    strT=['$\quad\quad$ $J^*(\mathbf{x})$ = ',sprintf('%10.3e',(min(minRes)))];
    strT={strT,['$J^*(\mathbf{x})-J^*_T$ = ',sprintf('%10.3e',min(minRes)-knownOptim)]};
    
    strT=regexprep(strT,'\ ','\\space');
    text(xT,yT,strT, 'interpreter','latex','HorizontalAlignment','right');
    
    
    % Plot 2
    axh=subplot(1,2,2,'ticklabelinterpreter','latex');
    
    switch dirOptim
        case 'min'
            errMeasure=-(knownOptim-minRes);
        case 'max'
            errMeasure=knownOptim-minRes;
    end
    
  
    
    
    popNominal=zeros([ceil(nIter/2),numel(optimstruct(1).population(1).fill)]);
    objNominal=zeros([ceil(nIter/2),1]);
    kk=0;
    for ii=1:2:nIter
        kk=kk+1;
        popNominal(kk,:)=optimstruct(ii).population(1).fill;
        objNominal(kk)=optimstruct(ii).population(1).objective;
    end
    
    
    iterVec=((1:2:nIter)+1)/2;
    varVec=1:length(popNominal(1,:));
    
    [iterGrid,varGrid]=meshgrid(iterVec,varVec);
    [objGrid,~]=meshgrid(objNominal,varVec);
    
    surf(iterGrid,varGrid,popNominal',objGrid)
    
    colorbar;
    
    xlabel('Iteration', 'interpreter','latex','fontsize',12)
    ylabel('Variable Index', 'interpreter','latex','fontsize',12)
    zlabel('Variable Value', 'interpreter','latex','fontsize',12)
    set(axh,'ticklabelinterpreter','latex')
    
    %% figure 2
    
    
    
    h(2)=figure('Name','Design Variable Evolution','Position',[20 100 1000 600]);
    
    axh=subplot(2,2,1,'ticklabelinterpreter','latex');
    s(1)=surf(iterGrid,varGrid,meanFill');
    caxis([0 1])
    
    axc(1)=colorbar;
    view(0 ,90)
    xlabel('Iteration', 'interpreter','latex','fontsize',12)
    ylabel('Variable Index', 'interpreter','latex','fontsize',12)
    
    
    
    cLabel={'Mean','Standard Deviation','Maximum','Minimum'};
    for ii=1:length(s)
        s(ii).EdgeColor='none';
        axc(ii).TickLabelInterpreter='latex';
        axc(ii).Label.String=cLabel{ii};
        axc(ii).Label.FontSize=14;
    end
    
    
end


function [iterRes,nIter,nVar]=BuildIterRes(optimstruct,defaultVal)
    
    nVar=0;
    for ii=1:length(optimstruct)
        nVar=max([length(optimstruct(ii).population),nVar]);
    end
    nIter=length(optimstruct);
    iterRes=zeros([nIter,nVar])+sign(defaultVal)*abs(defaultVal^2);
    
end
