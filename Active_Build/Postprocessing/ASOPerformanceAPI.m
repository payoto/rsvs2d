function [ASOstruct,h]=ASOPerformanceAPI(optIn, ASOiters,varargin)%dirSave,nameRun,figList
    % prepare for different cases
    runASOextract=1;
    
    [dirSave,nameRun,figList,splitCase]=HandleVarargin(varargin);
    isParam=false;
    isOptimStruct=false;
    if isstruct(optIn)
        try
            optimstruct=optIn.optimstruct;
            paramoptim=optIn.paramoptim;
            isOptimStruct=true;
        catch
            runASOextract=0;
            ASOstruct=optIn;
        end
    elseif ischar(optIn)
        try
            optPath=FindDir(optIn,'OptimRes',0);
            paramPath=FindDir(optIn,'FinalParam',0);
            optimstruct=load(optPath{1});
            optimstruct=optimstruct.optimstruct;
            paramoptim=load(paramPath{1});
            paramoptim=paramoptim.paramoptim;
            isParam=true;
            isOptimStruct=true;
        catch
            [optimstruct]=RecoverProfileLocations(optIn,ASOiters);
        end
        if isempty(dirSave)
            
            dirSave=optIn;
        end
        if isempty(nameRun)
            nameRun=regexprep(optIn,'^.*Dir[0-9,\-,_]*T[0-9]*_','');
        end
    elseif iscell(optIn)
        
        paramPath=FindDir(optIn{1},'FinalParam',0);
        
        paramoptim=load(paramPath{1});
        paramoptim=paramoptim.paramoptim;
        isParam=true;
        for ii=1:numel(optIn)
            optPath=FindDir(optIn{ii},'OptimRes',0);
            loadDat=load(optPath{1});
            if exist('optimstruct','var')
                optimstruct=[optimstruct,loadDat.optimstruct(ASOiters)];
                
            else
                optimstruct=loadDat.optimstruct;
            end
            
        end
        ASOiters=1:numel(optimstruct);
        isOptimStruct=true;
    else
        error('unknown input argument type')
    end
    % Do stuff
    knownOptim=0;
    
    kk=1;
    ll=0;
    if runASOextract
        if isParam
            varExtract={'direction','knownOptim','objectiveName'};
            [direction,knownOptim,objectiveName]=ExtractVariables(varExtract,paramoptim);
        end
        % ---------
        % This will cause issues when using many optimstructs stuck
        % together
        if isOptimStruct
            rmPop=false(size(optimstruct));
            for ii=1:numel(optimstruct)
                rmPop(ii)=isempty(optimstruct(ii).population(1).objective);
            end
            optimstruct=optimstruct(~rmPop);
        end
        % ---------
        for ii=ASOiters
            for jj=1:numel(optimstruct(ii).population)
                try
                    ASOstruct(kk)=ASOInterface(optimstruct(ii).population(jj).location);
                    ASOstruct(kk).DEIter=ii-1;
                    ASOstruct(kk).location=optimstruct(ii).population(jj).location;
                    kk=kk+1;
                catch MEid
                    MEid.getReport
                    disp([int2str(ii),':',int2str(jj),' failed'])
                    ll=ll+1;
                end
                
            end
        end
        disp([int2str(ll),' failures to ',int2str(kk),' success'])
        if isOptimStruct
            h=OptimHistory(0,optimstruct,knownOptim,1000,'min');
            ax=findobj(h(1),'type','axes');
        else
            h=figure('Name','ConvHistory');
            
            ax(2)=axes;
            hold on
        end
        
        if ~isempty(dirSave)
            
            save([dirSave,filesep,'ASOperformance',nameRun,'.mat'],'ASOstruct');
            assignin('base','ASOstruct',ASOstruct);
        end
    else
        h=figure;
        ax(2)=axes;
        hold on
    end
    
    [ASOstruct]=PlotASOPerformance_DataPreProc(ASOstruct);
    [h2]=PlotASOPerformance(ASOstruct,ax(2),splitCase);
    h=[h,h2];
    close(h(figList))
    % Save Data
    if ~isempty(dirSave)
        
        save([dirSave,filesep,'ASOperformance',nameRun,'.mat'],'ASOstruct');
        for ii=1:numel(h)
            h(ii).Name= [h(ii).Name,' - ',nameRun];
        end
        hgsave(h,[dirSave,filesep,'ASOperformance',nameRun,'.fig'])
        
    end
    
end

function [optimstruct]=RecoverProfileLocations(pathStr,iternum)
    
    for ii=iternum
        iterPaths=FindDirRegex(pathStr,['iteration_',int2str(ii),'$'],1);
        profPaths=FindDirRegex(iterPaths{1},'profile_',1);
        
        [optimstruct(ii).population(1:numel(profPaths)).location]...
            =deal(profPaths{:});
        
    end
end
    
function [dirSave,nameRun,figList,splitCase]=HandleVarargin(cellArgin)
    
    dirSave='';
    nameRun='';
    figList=[];
    splitCase='errorVecMode';
    
    for ii=1:2:numel(cellArgin)
        eval([cellArgin{ii},'=cellArgin{ii+1};']);
    end
    
end

function [ASOstruct]=ASOInterface(pathToASO)
    expectFields={'DEIter','majorIt','obj','opt','geomStepMag','eDV', 'loops',...
        'refLvl','geomErrMag','ASOdesVec','errX','errY','errNorm','errCNorm','errRaw',...
        'nSurfPoints','objFuncCalls','CD0','errorVecMode',...
        'location','nTopo','residual'};
    standardIn=cell(size(expectFields));
    [standardIn{:}]=deal(0);
    structBuild=[expectFields;standardIn];
    ASOstruct=struct(structBuild{:});
    
    [subdivLevel1, errorMagnitude, nDV,errX, errY, errNorm, errCNorm,errRaw] ...
        = ASO.Postproc.subdivData(pathToASO);
    [majorIt, objective, subdivLevel, geomStep,eDV, loops, opt, residual] = ...
        ASO.Postproc.iterationData(pathToASO);
    [CD0, nFunCall, nSurfPts,errorMode] = ASO.Postproc.ASOData(pathToASO);
    profiles = ASO.Postproc.profileData(pathToASO);
    
    ASOstruct.DEIter=1;
    ASOstruct.majorIt=majorIt;
    ASOstruct.obj=objective;
    ASOstruct.refLvl=(subdivLevel);
    ASOstruct.geomErrMag=errorMagnitude(FindObjNum([],unique(subdivLevel),(subdivLevel1)));
    ASOstruct.ASOdesVec=nDV(FindObjNum([],unique(subdivLevel),(subdivLevel1)))';
    ASOstruct.nSurfPoints=nSurfPts;
    ASOstruct.geomStepMag=geomStep;
    ASOstruct.objFuncCalls=nFunCall;
    ASOstruct.CD0=CD0;
    ASOstruct.errorVecMode=errorMode;
    ASOstruct.nTopo=numel(profiles);
    ASOstruct.errX =errX ;
    ASOstruct.errY=errY;
    ASOstruct.errNorm=errNorm;
    ASOstruct.errCNorm=errCNorm;
    ASOstruct.eDV=eDV;
    ASOstruct.loops=loops;
    ASOstruct.errRaw=errRaw;
    ASOstruct.opt=opt;
    ASOstruct.residual=residual;
    
    if isempty(majorIt)
        error('No ASO was run')
    end
    
    if ~isfinite(objective(1))
        error('ASO failed')
    end
    
end

function [h,ax]=PlotASOPerformance(ASOstruct,axDeOpt,splitCase,axOther)
    isaxDef=false;
    if nargin==4
        isaxDef=true;
        ax=axOther;
    end
    fieldsASO=fieldnames(ASOstruct);
    expectFields={'DEIter','majorIt','obj','refLvl','geomErrMag','ASOdesVec',...
        'nSurfPoints','geomStepMag','objFuncCalls','CD0'};
    %              scalar    vec       vec    vec     vec(numLvl)  vec(numLvl)    scalar     vec            scalar     scalar
    
    c=get(axDeOpt,'colororder');
    fColor=[ 0 0 0];
    ASOstructAll=ASOstruct;
    
    switch splitCase
        case 'errorVecMode'
            errVecModes={ASOstruct.errorVecMode};
            lErrVec=unique(errVecModes);
            funcLogTest=@(in) ~cellfun(@isempty,regexp(errVecModes,in));
        case 'subDivLevel'
            for ii=1:numel(ASOstruct)
                errVecModes{ii}=int2str(unique(ASOstruct(ii).refLvl)');
            end
            
            lErrVec=unique(errVecModes);
            funcLogTest=@(in) ~cellfun(@isempty,regexp(errVecModes,in));
        case 'RunName'
            errVecModes=regexprep(regexprep({ASOstruct.location}...
                ,'^.*Dir_[0-9-]*T[0-9]{6}_',''),'/.*$','');
            lErrVec=unique(errVecModes);
            funcLogTest=@(in) ~cellfun(@isempty,regexp(errVecModes,in));
        otherwise
            error('Unknown split case')
            
    end
    h2=figure('Name','Convergence History');
    ax2=axes(h2);
    hold on
    for iii=1:numel(lErrVec)
        ASOstruct=ASOstructAll(funcLogTest(lErrVec{iii}));
        
        if numel(lErrVec)>1
            fColor=c(mod(iii-1,size(c,1))+1,:);
        end
        % Iteration plots overlaid on DEresults
        for ii=1:numel(ASOstruct)
            plot(axDeOpt,[ASOstruct(ii).DEIter,ASOstruct(ii).DEIter+1],...
                [ASOstruct(ii).CD0,ASOstruct(ii).obj(end)],'--','color',fColor);

        end
        
        % Data rearrangement
        [summaryStruct(iii)]=PlotASOPerformance_DataExtraction(ASOstruct,lErrVec{iii});
        
        for ii=1:numel(ASOstruct)
            nums=cellfun(@str2double,regexp(regexprep(ASOstruct(ii).location,...
                '^.*profile_',''),'_','split'));
            plot3(ax2,ASOstruct(ii).majorIt+ASOstruct(ii).DEIter,...
                ones([1 numel(ASOstruct(ii).obj)])*nums(end),ASOstruct(ii).obj,'.-','color',fColor);
           
        end
    end
    view(ax2,0,0);
    % Generate figures
    if ~isaxDef
        figNames={'ASO Performance','ASO Performance Normalised',...
            'ASO desvar to surface','Func Maj','Scatter','Objective',...
            'Step Data','Box Plots'};
        figNamesSingle={};
        tfunc=@(x)numel(unique(x));
        runExtraFig=any(cellfun(tfunc,{summaryStruct.subDivLevel}));
        if runExtraFig
            
            kk=numel(figNames);
            kkStart=4*kk;
            xStr={'subDivLevel','nDesVarperBody','nonBasisDesignStart','nonBasisDesignFinal',...
                {'errorMeasure','errNormDiffdat','norm'},...
                {'errorMeasure','errNormdat','norm'},...
                {'errorMeasure','errCNormdat','norm'}};
            for jj=1:numel(xStr)
                kk=numel(figNames);
                if iscell(xStr{jj})
                    figNames(kk+(1:2))={['ASO Performance (',xStr{jj}{:},')'],...
                        ['ASO Performance Normalised (',xStr{jj}{:},')']};
                else
                    figNames(kk+(1:2))={['ASO Performance (',xStr{jj},')'],...
                        ['ASO Performance Normalised (',xStr{jj},')']};
                end
            end
            
        end
        
        kk=0;
        for ii=1:numel(figNames)
            h(ii)=figure('Name',figNames{ii});
            
            for jj=1:4
                kk=kk+1;
                ax(kk)=subplot(2,2,jj);
                hold on
            end
        end
        for ii=1:numel(figNamesSingle)
            h2(ii+1)=figure('Name',figNamesSingle{ii});
            axSingle(ii)=axes;
            hold on
            
        end
    end
    
    % Plot Data
    fColor=[ 0 0 0];
    fMarker='.:';
    markList='.+';
    lErrVec=regexprep(lErrVec,'_',' ');
    
    PlotASOPerformance_BoxPlot(ax(29),summaryStruct,lErrVec,...
        'changeObjRat','Obj Value change (% of orig objective)',0)
    PlotASOPerformance_BoxPlot(ax(30),summaryStruct,lErrVec,...
        'changeObj','Obj Value change',0)
    PlotASOPerformance_BoxPlot(ax(31),summaryStruct,lErrVec,...
        'bestObj','Best Objective',0.2)
    PlotASOPerformance_BoxPlot(ax(32),summaryStruct,lErrVec,...
        'geomStepMag','Geometric step',0)
    for iii=1:numel(summaryStruct)
        if numel(lErrVec)>1
            fColor=c(mod(iii-1,size(c,1))+1,:);
            if numel(lErrVec)>size(c,1)
                fMarker=markList(mod(floor((iii-1)/size(c,1)),size(markList,2))+1);
            end
        end
        PlotASOPerformance_fig1(ax(1:4),summaryStruct(iii),lErrVec{iii},fMarker,fColor)
        PlotASOPerformance_fig2(ax(5:8),summaryStruct(iii),lErrVec{iii},fMarker,fColor)
        PlotASOPerformance_fig3(ax(9:12),summaryStruct(iii),lErrVec{iii},fMarker,fColor)
        PlotASOPerformance_fig4(ax(13:16),summaryStruct(iii),lErrVec{iii},fMarker,fColor)
        PlotASOPerformance_fig5(ax(17:20),summaryStruct(iii),lErrVec{iii},fMarker,fColor)
        PlotASOPerformance_fig6(ax(21:24),summaryStruct(iii),lErrVec{iii},fMarker,fColor)
        PlotASOPerformance_fig7(ax(25:28),summaryStruct(iii),lErrVec{iii},fMarker,fColor)
        % Variable subdiv plots
        kk=kkStart;
        for jj=1:numel(xStr)
            
            if iscell(xStr{jj})
                x=summaryStruct(iii).(xStr{jj}{1});
                for ll=2:numel(xStr{jj})
                    x=x.(xStr{jj}{ll});
                end
                xName=[xStr{jj}{:}];
            else
                x=summaryStruct(iii).(xStr{jj});
                xName=xStr{jj};
            end
            if runExtraFig
                PlotASOPerformance_asoperffigx(ax(kk+(1:4)),x,xName,...
                    summaryStruct(iii),lErrVec{iii},fMarker,fColor)
                PlotASOPerformance_asoperffignormx(ax(kk+(5:8)),x,xName,...
                    summaryStruct(iii),lErrVec{iii},fMarker,fColor)
                kk=kk+8;
            end
        end
        
    end
    
    for ii=[2 3 5 6 7  9 11 16 21 23 ...
            kkStart+(9:16) ...
            kkStart+2*8+(1:8)  kkStart+3*8+(1:8) kkStart+5*8+(1:8) kkStart+6*8+(1:8)]
        ax(ii).XScale='log';
    end
    % err mag plots
    AverageAllLines(ax(2),[1e-3 1e-2 1e-1 1 10])
    AverageAllLines(ax(3),[1e-3 1e-2 1e-1 1 10])
    AverageAllLines(ax(6),[1e-3 1e-2 1e-1 1 10])
    AverageAllLines(ax(7),[1e-3 1e-2 1e-1 1 10])
    AverageAllLines(ax(9),{'lin',4})
    AddNumbersInRange(ax(3),0)
    AddNumbersInRange(ax(17),0.9,2)
    AddNumbersInRange(ax(18),0.9,2)
    AddNumbersInRange(ax(19),0.9,2)
    AddNumbersInRange(ax(20),0.9,2)
    %AddNumbersInRange(ax(12),0)
    
    AverageAllLines(ax(12),{{'lin',4} , {'lin',4}})
    AverageAllLines(ax(4),{'lin',4})
    
    AverageAllLines(ax(kkStart+7),[1.5:4.5])
    AverageAllLines(ax(10),{750:250:2000 , 50:50:250})
    AverageAllLines(ax(13),{'lin',3})
    AverageAllLines(ax(14),{'lin',3})
    AverageAllLines(ax(15),{'lin',3})
    AverageAllLines(ax(16),{'log',3})
    AverageAllLines(ax(22),[1:10:101])
    AverageAllLines(ax(24),[1:10:101])
    AverageAllLines(ax(21),[0.03:0.02:0.11])
    AverageAllLines(ax(23),{'log',5})
    
    AverageAllLines(ax(28),[1.5:4.5])
    % sub div
    
    AverageAllLines(ax(kkStart+2),[1.5:4.5])
    AverageAllLines(ax(kkStart+3),[1.5:4.5])
    AverageAllLines(ax(kkStart+6),[1.5:4.5])
    AverageAllLines(ax(kkStart+7),[1.5:4.5])
    AverageAllLines(ax(kkStart+10),{'log',3})
    AverageAllLines(ax(kkStart+11),{'log',3})
    AverageAllLines(ax(kkStart+12),{'log',3})
    AverageAllLines(ax(kkStart+14),{'log',3})
    AverageAllLines(ax(kkStart+15),{'log',3})
    AverageAllLines(ax(kkStart+16),{'log',3})
    for ii=kkStart+4*8+1:kkStart+5*8
        AverageAllLines(ax(ii),{'lin',4})
    end
    for ii=[kkStart+2*8+1:kkStart+4*8,kkStart+5*8+1:numel(ax)]
        AverageAllLines(ax(ii),{'log',4})
    end
    legend(ax(1),findobj(ax(1),'type','line'))
    PlotSummaryStruct(summaryStruct)
    
    h=[h2,h];
end

function [ASOstruct]=PlotASOPerformance_DataPreProc(ASOstruct)
    
    normVec=@(x) sqrt(sum(x.^2,2));
    
    for ii=1:numel(ASOstruct)
        % Calculate a value for the design represented by the error vectors
        ASOstruct(ii).nonBasisDesign=zeros(size(ASOstruct(ii).obj));
        actRef=unique(ASOstruct(ii).refLvl);
        errX=[ASOstruct(ii).errX{actRef}];
        errY=[ASOstruct(ii).errY{actRef}];
        errMat=cell(size(actRef));
        for jj=1:numel(actRef)
            errMat{jj}=[errX(:,jj),errY(:,jj),errY(:,jj),errX(:,jj)];
        end
        for jj=1:numel(ASOstruct(ii).nonBasisDesign)
            errVec=[errMat{ASOstruct(ii).refLvl(jj)==...
                actRef}(:,[1 3])*ASOstruct(ii).eDV(jj,[1 3])', ...
                errMat{ASOstruct(ii).refLvl(jj)==...
                actRef}(:,[2 4])*ASOstruct(ii).eDV(jj,[2 4])'];
            ASOstruct(ii).nonBasisDesign(jj)=sqrt(sum(normVec(errVec).^2));
            ASOstruct(ii).nonBasisDesignX(jj)=sqrt(sum((errVec(:,1)).^2));
            ASOstruct(ii).nonBasisDesignY(jj)=sqrt(sum((errVec(:,2)).^2));
        end
       
        % Calculate Measures for the onesidedness of the error Vector
        for jj=1:numel(actRef)
            xy=normVec([ASOstruct(ii).errX{actRef(jj)},...
                ASOstruct(ii).errY{actRef(jj)}]);
            ASOstruct(ii).errXYdat(jj)=...
                GetErrorProperties(normVec([ASOstruct(ii).errX{actRef(jj)},...
                ASOstruct(ii).errY{actRef(jj)}]));
            ASOstruct(ii).errNormdat(jj)=...
               GetErrorProperties( ASOstruct(ii).errNorm{actRef(jj)});
            ASOstruct(ii).errCNormdat(jj)=...
                GetErrorProperties(ASOstruct(ii).errCNorm{actRef(jj)});
            ASOstruct(ii).errNormDiffdat(jj)=...
                GetErrorProperties(real(sqrt(xy.^2-...
                ASOstruct(ii).errNorm{actRef(jj)}.^2))/ASOstruct(ii).geomErrMag(jj));
            if any(imag(ASOstruct(ii).errNormDiffdat(jj).norm)>0)
                warning('is imaginary')
            end
        end
    end
    
end

function [datstruct]=GetErrorProperties(vec)
    
    datstruct.norm=sqrt(sum(vec.^2));
    
    datstruct.mean=mean(vec);
    datstruct.median=median(vec);
    datstruct.std=std(vec);
    
    datstruct.medianabs=median(abs(vec));
    datstruct.meanabs=mean(abs(vec));
    datstruct.stdabs=std(abs(vec));
    
    datstruct.delta=max(vec)-min(vec);
    datstruct.max=max(vec);
    datstruct.min=min(vec);
    datstruct.centre=(max(vec)+min(vec))/2;
    
    
end

function [summaryStruct]=PlotASOPerformance_DataExtraction(ASOstruct,lErrVec)
    
    fieldsFullExplore={'errXYdat','errNormdat','errCNormdat','errNormDiffdat'};
    if true
        ii=1;
        jj=1;
        d2=[];
        keepVar=[];
        d=[];
        d=who;
    end
    errMagVec=[];    subDivLevel=[];    changeObj=[];    changeObjRat=[];
    changeObjStep=[];    changeObjStepRat=[];
    nIter95Delta=[];    nIter95DeltaRat=[];    nDesVar=[];    nDesVarperBody=[];
    changeObjInitial=[];    changeObjInitialRat=[];    changeObjCD0=[];
    changeObjFinalCD0=[];    geomStepMag=[];    nSurfPoints=[];
    nFuncperMaj=[];   changeAtIter10=[]; changeAtIter10Rat=[]; bestObj=[];
    cd0=[]; profNum=[];nonBasisDesignFinal=[];nonBasisDesignStart=[];
    changeAtIter10LocalLvl=[];changeAtIter10LocalLvlRat=[];
    
    for ii=fieldsFullExplore
        for jj=fieldnames(ASOstruct(1).(ii{1}))'
            eval(['errorMeasure.',ii{1},'.',jj{1},'=[];']);
        end
    end
    
    if true
        d2=who;
        keepVar=true(size(d2));
        for ii=1:numel(d)
            keepVar=keepVar & cellfun(@(x)~strcmp(x,d{ii}),d2);
        end
        keepCell=[d2(keepVar);'errorVecMode']';
        standardIn=cell(size(keepCell));
        varsAct=[d2(keepVar)];
        
        structBuild=[keepCell;standardIn];
        summaryStruct=struct(structBuild{:});
    end
    
    deltaPercent=0.9;
    iter10=10;
    
    for ii=1:numel(ASOstruct)
        jjStart=numel(errMagVec);
        errMagVec=[errMagVec,ASOstruct(ii).geomErrMag'];
        nDesVar=[nDesVar,ASOstruct(ii).ASOdesVec];
        nDesVarperBody=[nDesVarperBody,ASOstruct(ii).ASOdesVec/ASOstruct(ii).nTopo];
        subDivLevel=[subDivLevel,unique(ASOstruct(ii).refLvl)'];
        cd0=[cd0,unique(ASOstruct(ii).CD0)*ones(size(unique(ASOstruct(ii).refLvl)))'];
        nums=cellfun(@str2double,regexp(regexprep(ASOstruct(ii).location,...
            '^.*profile_',''),'_','split'));
        profNum=[profNum,nums(end)*ones(size(unique(ASOstruct(ii).refLvl)))'];
        
        % Bulk variables
        for ii1=fieldsFullExplore
            for jj1=fieldnames(ASOstruct(ii).(ii1{1}))'
                varName=['errorMeasure.',ii1{1},'.',jj1{1}];
                eval([varName,'=[',varName,',ASOstruct(ii).(ii1{1}).(jj1{1})];']);
            end
        end
        
        
        if numel(errMagVec)~=numel(subDivLevel)
            error('subDivLevel size does not match number of errMagnitudes provided')
        end
        for jj=1:(numel(errMagVec)-jjStart)
            currList=find(ASOstruct(ii).refLvl==subDivLevel(jj+jjStart));
            if numel(currList)==1
                currList(2)=currList(1);
            end
            bestObj(jjStart+jj)=min(ASOstruct(ii).obj(currList));
            changeObj(jjStart+jj)=min(ASOstruct(ii).obj(currList))...
                -ASOstruct(ii).obj(1);
            changeObjRat(jjStart+jj)=changeObj(jjStart+jj)/ASOstruct(ii).obj(1); %#ok<*AGROW>
            
            changeObjStep(jjStart+jj)=min(ASOstruct(ii).obj(currList))...
                -ASOstruct(ii).obj(currList(1));
            changeObjStepRat(jjStart+jj)=changeObjStep(jjStart+jj)/ASOstruct(ii).obj(currList(1));
            
            changeObjInitial(jjStart+jj)=(ASOstruct(ii).obj(currList(2))...
                -ASOstruct(ii).obj(currList(1)));
            changeObjInitialRat(jjStart+jj)=(ASOstruct(ii).obj(currList(2))...
                -ASOstruct(ii).obj(currList(1)))/(ASOstruct(ii).obj(currList(end))...
                -ASOstruct(ii).obj(currList(1)));
            
            if changeObj(jjStart+jj)<0
                nIter95Delta(jjStart+jj)=min(find((ASOstruct(ii).obj(1)+...
                    changeObj(jjStart+jj)*deltaPercent)>=ASOstruct(ii).obj(currList)));
            else
                nIter95Delta(jjStart+jj)=0;
            end
            nIter95DeltaRat(jjStart+jj)=nIter95Delta(jjStart+jj)/numel(currList);
            nSurfPoints(jj+jjStart)=ASOstruct(ii).nSurfPoints;
            
            
            changeObjCD0(jj+jjStart)=ASOstruct(ii).obj(1)-ASOstruct(ii).CD0;
            changeObjFinalCD0(jj+jjStart)=ASOstruct(ii).obj(end)-ASOstruct(ii).CD0;
            geomStepMag(jj+jjStart)=ASOstruct(ii).geomStepMag(currList(end));
            nFuncperMaj(jj+jjStart)=ASOstruct(ii).objFuncCalls/numel(currList);
            changeAtIter10(jj+jjStart)=(ASOstruct(ii).obj(currList(...
                min(iter10,numel(currList))))-ASOstruct(ii).obj(1));
            changeAtIter10Rat(jj+jjStart)=changeAtIter10(jj+jjStart)/changeObj(jjStart+jj);
            
            changeAtIter10LocalLvl(jj+jjStart)=(ASOstruct(ii).obj(currList(...
                min(iter10,numel(currList))))-ASOstruct(ii).obj(currList(...
                min(1,numel(currList)))));
            changeAtIter10LocalLvlRat(jj+jjStart)=changeAtIter10LocalLvl(jj+jjStart)/changeObj(jjStart+jj);
            nonBasisDesignFinal(jj+jjStart)=ASOstruct(ii).nonBasisDesign(currList(end));
            nonBasisDesignStart(jj+jjStart)=ASOstruct(ii).nonBasisDesign(currList(1));
        end
        % Add nans in between subdivision runs to see them
        for ll=1:numel(varsAct)
            if ~strcmp(varsAct{ll},'errorMeasure')
                eval([varsAct{ll},'=[',varsAct{ll},',nan];']);
            else
                for ii1=fieldsFullExplore
                    for jj1=fieldnames(ASOstruct(ii).(ii1{1}))'
                        varName=['errorMeasure.',ii1{1},'.',jj1{1}];
                        eval([varName,'=[',varName,',nan];']);
                    end
                end
            end
        end
    end
    fieldsSum=fieldnames(summaryStruct);
    for ii=1:numel(fieldsSum)-1
        summaryStruct.(fieldsSum{ii})=eval(fieldsSum{ii});
    end
    summaryStruct.(fieldsSum{end})=lErrVec;
    
end

function []=PlotASOPerformance_BoxPlot(ax,summarrystruct,caseName,fieldName,...
        str,replaceNan)
    if ~exist('str','var');str=fieldName;end
    if ~exist('replaceNan','var');replaceNan=nan;end
    
    nDat=cellfun(@numel,{summarrystruct.(fieldName)});
    for ii=1:numel(summarrystruct)
        for jj=1:nDat(ii)-1
            if ~isnan(summarrystruct(ii).(fieldName)(jj+1))
                summarrystruct(ii).(fieldName)(jj)=nan;
            end
        end
        summarrystruct(ii).(fieldName)(isnan(summarrystruct(ii).(fieldName)))=[];
    end
    
    nDat=cellfun(@numel,{summarrystruct.(fieldName)});
    boxDat=nan([max(nDat),numel(summarrystruct)]);
    boxDat(isnan(boxDat))=replaceNan;
    for ii=1:numel(summarrystruct)
        boxDat(1:nDat(ii),ii)=summarrystruct(ii).(fieldName);
    end
    
    
    boxplot(ax,boxDat,caseName)
    
    ax.YLabel.String=str; 
end

function []=PlotASOPerformance_BoxPlotLog(ax,summarrystruct,caseName,fieldName,str)
    if ~exist('str','var');str=fieldName;end
    
    nDat=cellfun(@numel,{summarrystruct.(fieldName)});
    for ii=1:numel(summarrystruct)
        for jj=1:nDat(ii)-1
            if ~isnan(summarrystruct(ii).(fieldName)(jj+1))
                summarrystruct(ii).(fieldName)(jj)=nan;
            end
        end
        summarrystruct(ii).(fieldName)(isnan(summarrystruct(ii).(fieldName)))=[];
    end
    
    nDat=cellfun(@numel,{summarrystruct.(fieldName)});
    boxDat=nan([max(nDat),numel(summarrystruct)]);
    
    for ii=1:numel(summarrystruct)
        boxDat(1:nDat(ii),ii)=summarrystruct(ii).(fieldName);
    end
    
    
    boxplot(ax,log10(boxDat),caseName)
    
    ax.YLabel.String=str; 
end

function []=PlotASOPerformance_fig1(ax,summarrystruct,caseName,fMarker,fColor)
    
    
    % Error Magnitude (and later subdivision level) vs
    % 1) geometric error Magnitude (cloud)
    plot(ax(1),summarrystruct.subDivLevel,summarrystruct.errMagVec,fMarker,'color',...
        fColor,'DisplayName',caseName)
    % 2) Change in objective thanks to ASO (cloud)
    plot(ax(2),summarrystruct.errMagVec,summarrystruct.changeObj,fMarker,'color',fColor)
    % 3) Initial Objective Delta (should be 0 except for error = 'none')
    plot(ax(3),summarrystruct.errMagVec,summarrystruct.changeObjInitial,fMarker,'color',fColor)
    % 4) Change of objective other the first iteration (normalised by total Delta)
    plot(ax(4),summarrystruct.nDesVar,summarrystruct.nIter95Delta,fMarker,'color',fColor)
    
    
    ax(1).XLabel.String='subdivision level';
    ax(1).YLabel.String='subdiv error mag';
    ax(2).XLabel.String='subdiv error mag';
    ax(2).YLabel.String='Obj Value change';
    ax(3).XLabel.String='subdiv error mag';
    ax(3).YLabel.String='Obj Value change over majorIt 1';
    ax(4).XLabel.String='# Design Variables';
    ax(4).YLabel.String='nIter to 90% obj change';
    
end

function []=PlotASOPerformance_fig2(ax,summarrystruct,caseName,fMarker,fColor)
    
    % Error Magnitude (and later subdivision level) vs
    % 1) geometric error Magnitude (cloud)
    plot(ax(1),summarrystruct.errMagVec,summarrystruct.nDesVar,fMarker,'color',fColor)
    % 2) Change in objective thanks to ASO (cloud)
    plot(ax(2),summarrystruct.errMagVec,summarrystruct.changeObjRat,fMarker,'color',fColor)
    % 3) Initial Objective Delta (should be 0 except for error = 'none')
    plot(ax(3),summarrystruct.errMagVec,summarrystruct.changeObjInitialRat,fMarker,'color',fColor)
    % 4) Change of objective other the first iteration (normalised by total Delta)
    plot(ax(4),summarrystruct.nDesVar,summarrystruct.nIter95DeltaRat,fMarker,'color',fColor)
    
    
    
    ax(1).XLabel.String='subdiv error mag';
    ax(1).YLabel.String='# Design Variables';
    ax(2).XLabel.String='subdiv error mag';
    ax(2).YLabel.String='Obj Value change (% of orig objective)';
    ax(3).XLabel.String='subdiv error mag';
    ax(3).YLabel.String='Obj Value change over majorIt 1 (% total change)';
    ax(4).XLabel.String='# Design Variables';
    ax(4).YLabel.String='nIter to 90% obj change (% total #iter)';
    
end

function []=PlotASOPerformance_fig3(ax,summarrystruct,caseName,fMarker,fColor)
    
    %plot(ax(1),summarrystruct.nSurfPoints,summarrystruct.nDesVar,fMarker,'color',fColor)
    plot(ax(1),summarrystruct.nDesVar,summarrystruct.geomStepMag,fMarker,'color',fColor)
    % Magnitude of the geometric step vs number of surf points and # of design variables
    plot3(ax(2),summarrystruct.nSurfPoints,summarrystruct.nDesVar,...
        summarrystruct.geomStepMag,fMarker,'color',fColor)
    semilogy(ax(3),summarrystruct.errMagVec,abs(summarrystruct.changeObjCD0),fMarker,'color',fColor)
    ax(3).YScale='log';
    %plot(ax(3),summarrystruct.changeAtIter10Rat,summarrystruct.changeObj,fMarker,'color',fColor)
    
    plot3(ax(4),summarrystruct.nDesVarperBody,summarrystruct.changeObj,...
        summarrystruct.nIter95Delta,fMarker,'color',fColor)
    %plot(ax(4),summarrystruct.errMagVec,summarrystruct.changeObjFinalCD0,fMarker,'color',fColor)
    % Add lines and numbers indicating overUnder
    
    
    %     ax(1).XLabel.String='# surf points';
    %     ax(1).YLabel.String='# design variable';
    ax(1).XLabel.String='# desVar per body';
    ax(1).YLabel.String='Geometric step';
    ax(2).XLabel.String='# surf points';
    ax(2).YLabel.String='# design variable';
    ax(2).ZLabel.String='Size of Geometric Step';
    ax(3).XLabel.String='subdiv error mag';
    ax(3).YLabel.String='CD0-obj(1)';
    ax(4).XLabel.String='Des var per body';
    ax(4).YLabel.String='change in obj';
    ax(4).ZLabel.String='nIter 90%';
end

function []=PlotASOPerformance_fig4(ax,summarrystruct,caseName,fMarker,fColor)
    
    plot(ax(1),summarrystruct.subDivLevel,summarrystruct.nFuncperMaj,fMarker,'color',fColor)
    % Magnitude of the geometric step vs number of surf points and # of design variables
    plot(ax(2),summarrystruct.nDesVar,summarrystruct.nFuncperMaj,fMarker,'color',fColor)
    plot(ax(3),summarrystruct.nDesVarperBody,summarrystruct.nFuncperMaj,fMarker,'color',fColor)
    plot(ax(4),summarrystruct.errMagVec,summarrystruct.nFuncperMaj,fMarker,'color',fColor)
    % Add lines and numbers indicating overUnder
    
    
    ax(1).XLabel.String='# surf points';
    ax(1).YLabel.String='obj func/ majIter';
    ax(2).XLabel.String='# design variables';
    ax(2).YLabel.String='obj func/ majIter';
    ax(3).XLabel.String='# des var per topo';
    ax(3).YLabel.String='obj func/ majIter';
    ax(4).XLabel.String='subdiv error mag';
    ax(4).YLabel.String='obj func/ majIter';
end

function []=PlotASOPerformance_fig5(ax,summarrystruct,caseName,fMarker,fColor)
    
    s=scatter(ax(1),summarrystruct.changeAtIter10Rat,summarrystruct.changeObj,9,...
        summarrystruct.subDivLevel);
    ax(1).Color=[0.7 0.7 0.7];
    s.MarkerEdgeColor=fColor;
    s.MarkerFaceColor='flat';
    % Magnitude of the geometric step vs number of surf points and # of design variables
    s=scatter(ax(2),summarrystruct.changeAtIter10Rat,summarrystruct.changeObj,9,...
        log10(summarrystruct.nDesVarperBody));
    ax(2).Color=[0.7 0.7 0.7];
    s.MarkerEdgeColor=fColor;
    s.MarkerFaceColor='flat';
    %plot(ax(3),summarrystruct.errMagVec,summarrystruct.changeObjCD0,fMarker,'color',fColor)
    %plot(ax(3),summarrystruct.changeAtIter10Rat,summarrystruct.changeObj,fMarker,'color',fColor)
    s=scatter(ax(3),summarrystruct.changeAtIter10Rat,summarrystruct.changeObjRat,9,...
        (summarrystruct.subDivLevel));
    ax(3).Color=[0.6 0.6 0.6];
    s.MarkerEdgeColor=fColor;
    s.MarkerFaceColor='flat';
    s=scatter(ax(4),summarrystruct.changeAtIter10Rat,summarrystruct.changeObjRat,9,...
        log10(summarrystruct.nDesVarperBody));
    ax(4).Color=[0.6 0.6 0.6];
    s.MarkerEdgeColor=fColor;
    s.MarkerFaceColor='flat';
    %plot(ax(4),summarrystruct.errMagVec,summarrystruct.changeObjFinalCD0,fMarker,'color',fColor)
    % Add lines and numbers indicating overUnder
    
    
    [ax.XLim]=deal([0 1]);
    ax(1).XLabel.String='Ratio change at iter 10/total change';
    ax(1).YLabel.String='Total change';
    ax(1).Title.String='subDivLevel';
    ax(2).XLabel.String='Ratio change at iter 10/total change';
    ax(2).YLabel.String='Total change';
    ax(2).Title.String='nDesVarperBody';
    ax(3).XLabel.String='Ratio change at iter 10/total change';
    ax(3).YLabel.String='Total change/obj val';
    ax(3).Title.String='subDivLevel';
    ax(4).XLabel.String='Ratio change at iter 10/total change';
    ax(4).YLabel.String='Total change/obj val';
    ax(4).Title.String='nDesVarperBody';
end

function []=PlotASOPerformance_fig6(ax,summarrystruct,caseName,fMarker,fColor)
    
    
    % Rain plots
    %{
    l=plot(ax(1),[summarrystruct.cd0;summarrystruct.cd0*1.05],...
        [summarrystruct.cd0;summarrystruct.bestObj]);
    [l.Marker]=deal(fMarker);
    [l.Color]=deal(fColor);
    ax(1).YScale='log';
    l=plot(ax(2),[summarrystruct.profNum;summarrystruct.profNum+0.5],...
        [summarrystruct.cd0;summarrystruct.bestObj]);
    [l.Marker]=deal(fMarker);
    [l.Color]=deal(fColor);
    ax(2).YScale='log';
    %}
    ReorganiseSubPlots(ax,[1 4],[0.075 0.01 0.1 0.01],[0 0.075],[9 11])
    
    notNan=find(~isnan(summarrystruct.profNum));
    x2=summarrystruct.errMagVec(notNan);
    [~,i]=sort(x2);
    l=plot(ax(1),[x2(i)],[summarrystruct.changeObjRat(notNan(i))]);
    [l.Marker]=deal(fMarker(1));
    [l.Color]=deal(fColor);
    %ax(1).YScale='log';
    x=1:100;
    y=zeros(size(x));
    y([summarrystruct.profNum(notNan)])=[summarrystruct.changeObjRat(notNan)];
    l=plot(ax(2),x,y);
    [l.Marker]=deal(fMarker(1));
    [l.Color]=deal(fColor);
    %ax(2).YScale='log';
    % Line plots
    l=plot(ax(3),[x2(i);x2(i)]',...
        [summarrystruct.cd0(i);summarrystruct.bestObj(i)]');
    [l.Marker]=deal(fMarker(1));
    [l.Color]=deal(fColor);
    ax(3).YScale='log';
    
    y=zeros(size(x))*nan;
    y([summarrystruct.profNum(notNan)])=[summarrystruct.bestObj(notNan)];
    l=plot(ax(4),x,y);
    [l.Marker]=deal(fMarker(1));
    [l.Color]=deal(fColor);
    ax(4).YScale='log';
    
    
    ax(1).XLabel.String='Error Mag Vec ';
    ax(1).YLabel.String='Change of Objective';
    ax(2).XLabel.String='Profile number';
    ax(2).YLabel.String='Change of Objective';
    ax(3).XLabel.String='Error Mag Vec ';
    ax(3).YLabel.String='Best objective';
    ax(4).XLabel.String='Profile number';
    ax(4).YLabel.String='Best objective';
end

function []=PlotASOPerformance_fig7(ax,summarrystruct,caseName,fMarker,fColor)
    
    
    % Rain plots
    %{
    l=plot(ax(1),[summarrystruct.cd0;summarrystruct.cd0*1.05],...
        [summarrystruct.cd0;summarrystruct.bestObj]);
    [l.Marker]=deal(fMarker);
    [l.Color]=deal(fColor);
    ax(1).YScale='log';
    l=plot(ax(2),[summarrystruct.profNum;summarrystruct.profNum+0.5],...
        [summarrystruct.cd0;summarrystruct.bestObj]);
    [l.Marker]=deal(fMarker);
    [l.Color]=deal(fColor);
    ax(2).YScale='log';
    %}
    ReorganiseSubPlots(ax(1:2),[1 2],[0.075 0.01 0.51 0.01],[0 0.075],[9 11])
    
    notNan=find(~isnan(summarrystruct.profNum));
    x2=summarrystruct.errMagVec(notNan);
    [~,i]=sort(x2);
    x=1:100;
    y=zeros(size(x));
    y([summarrystruct.profNum(notNan)])=[summarrystruct.changeObj(notNan)];
    l=plot(ax(1),x,y);
    plot(ax(1),summarrystruct.profNum,summarrystruct.changeObj,fMarker,'color',fColor)
    [l.Marker]=deal(fMarker(1));
    [l.Color]=deal(fColor);
    %ax(1).YScale='log';
    x=1:100;
    y=zeros(size(x));
    y([summarrystruct.profNum(notNan)])=[summarrystruct.changeObjRat(notNan)];
    l=plot(ax(2),x,y);
    plot(ax(2),summarrystruct.profNum,summarrystruct.changeObjRat,fMarker,'color',fColor)
    [l.Marker]=deal(fMarker(1));
    [l.Color]=deal(fColor);
    %ax(2).YScale='log';
    % Line plots
    l=plot(ax(3),[x2(i);x2(i)]',...
        [summarrystruct.cd0(i);summarrystruct.bestObj(i)]');
    [l.Marker]=deal(fMarker(1));
    [l.Color]=deal(fColor);
    ax(3).YScale='log';
    
    y=zeros(size(x))*nan;
    plot(ax(4),summarrystruct.subDivLevel,summarrystruct.changeObjStepRat,fMarker,'color',fColor)
    ax(4).YScale='log';
    
    
    ax(1).XLabel.String='Error Mag Vec ';
    ax(1).YLabel.String='Change of Objective';
    ax(2).XLabel.String='Profile number';
    ax(2).YLabel.String='Change of Objective';
    ax(3).XLabel.String='Error Mag Vec ';
    ax(3).YLabel.String='Best objective';
    ax(4).XLabel.String='Subdivision Level';
    ax(4).YLabel.String='Change of objective per step';
end


function []=PlotASOPerformance_asoperffigx(ax,x,xStr,summarrystruct,caseName,fMarker,fColor)
    
    
    % Error Magnitude (and later subdivision level) vs
    % 1) geometric error Magnitude (cloud)
    plot(ax(1),x,summarrystruct.errMagVec,fMarker,'color',...
        fColor,'DisplayName',caseName)
    % 2) Change in objective thanks to ASO (cloud)
    plot(ax(2),x,summarrystruct.changeObj,fMarker,'color',fColor)
    % 3) Initial Objective Delta (should be 0 except for error = 'none')
    plot(ax(3),x,summarrystruct.changeObjInitial,fMarker,'color',fColor)
    % 4) Change of objective other the first iteration (normalised by total Delta)
    plot(ax(4),x,summarrystruct.nIter95Delta,fMarker,'color',fColor)
    
    for ii=1:4
        [ax(ii).XLabel.String]=(xStr);
    end
    ax(1).YLabel.String='subdiv error mag';
    ax(2).YLabel.String='Obj Value change';
    ax(3).YLabel.String='Obj Value change over majorIt 1';
    ax(4).YLabel.String='nIter to 90% obj change';
    
end

function []=PlotASOPerformance_asoperffignormx(ax,x,xStr,summarrystruct,caseName,fMarker,fColor)
    
    % Error Magnitude (and later subdivision level) vs
    % 1) geometric error Magnitude (cloud)
    plot(ax(1),x,summarrystruct.nDesVar,fMarker,'color',fColor)
    % 2) Change in objective thanks to ASO (cloud)
    plot(ax(2),x,summarrystruct.changeObjRat,fMarker,'color',fColor)
    % 3) Initial Objective Delta (should be 0 except for error = 'none')
    plot(ax(3),x,summarrystruct.changeObjInitialRat,fMarker,'color',fColor)
    % 4) Change of objective other the first iteration (normalised by total Delta)
    plot(ax(4),x,summarrystruct.nIter95DeltaRat,fMarker,'color',fColor)
    
    for ii=1:4
        [ax(ii).XLabel.String]=(xStr);
    end
    ax(1).YLabel.String='# Design Variables';
    ax(2).YLabel.String='Obj Value change (% of orig objective)';
    ax(3).YLabel.String='Obj Value change over majorIt 1 (% total change)';
    ax(4).YLabel.String='nIter to 90% obj change (% total #iter)';
    
end

function []=PlotGeomPerformance_asoperffig()
    
end

function []=PlotSummaryStruct(summaryStruct)
    
end


function []=AverageAllLines(ax,ranges)
    
    
    hold(ax,'on');
    isType=1;
    l=findobj(ax,'type','line');
    if isempty(l)
        l=findobj(ax,'type','scatter');
        isType=2;
    end
    for ii=1:numel(l)
        if isempty(l(ii).ZData)
            [laverage]=AverageLines(ax,l(ii).XData,l(ii).YData,ranges);
            if ~isempty(laverage)
                if isType==1
                    [laverage.Color]=deal(l(ii).Color);
                else
                    [laverage.Color]=deal(l(ii).MarkerEdgeColor);
                end
            end
        else
            [laverage]=AverageLines3D(ax,l(ii).XData,l(ii).YData,...
                l(ii).ZData,ranges);
            [laverage.EdgeColor]=deal(l(ii).Color);
        end
        if ~isempty(laverage)
            [laverage.LineStyle]=deal('-');
            laverage(1).Marker='*';
            laverage(2).Marker='d';
        end
    end
end

function [l]=AverageLines(ax,x,y,ranges)
    
    if iscell(ranges)
        switch ranges{1}
            case 'lin'
                ranges=linspace(min(x),max(x)*1.01,ranges{2}+1);
            case 'log'
                ranges=logspace(log10(min(x)),log10(max(x)*1.01),ranges{2}+1);
        end
        
    end
    
    hold(ax,'on');
    cellRanges=cell(numel(ranges)+1,1);
    ranges=[-inf,ranges,inf];
    xPlotMean=[];
    yPlotMean=[];
    xPlotMedian=[];
    yPlotMedian=[];
    for ii=1:numel(cellRanges)
        cellRanges{ii}=find(x<ranges(ii+1) & x>=ranges(ii));
        if numel(cellRanges{ii})>0
            xPlotMean=[xPlotMean,mean(x(cellRanges{ii}),'omitnan' )];
            yPlotMean=[yPlotMean,mean(y(cellRanges{ii}),'omitnan' )];
            xPlotMedian=[xPlotMedian,median(x(cellRanges{ii}),'omitnan' )];
            yPlotMedian=[yPlotMedian,median(y(cellRanges{ii}),'omitnan' )];
        end
    end
    
    l=plot(ax,xPlotMean,yPlotMean,xPlotMedian,yPlotMedian);
end

function [l]=AverageLines3D(ax,x,y,z,ranges)
    if iscell(ranges{1})
        switch ranges{1}{1}
            case 'lin'
                rangesX=linspace(min(x),max(x),ranges{1}{2}+1);
            case 'log'
                rangesX=logspace(log10(min(x)),log10(max(x)),ranges{1}{2}+1);
        end
        rangesX=rangesX(2:end-1);
        
    else
        rangesX=ranges{1};
        rangesY=ranges{2};
        
    end
    if iscell(ranges{2})
        switch ranges{2}{1}
            case 'lin'
                rangesY=linspace(min(y),max(y),ranges{2}{2}+1);
            case 'log'
                rangesY=logspace(log10(min(y)),log10(max(y)),ranges{2}{2}+1);
        end
        rangesY=rangesY(2:end-1);
        
    end
    
    hold(ax,'on');
    cellRanges=cell(numel(rangesX)+1,numel(rangesY)+1);
    rangesX=[-inf,rangesX,inf];
    rangesY=[-inf,rangesY,inf];
    
    xPlotMean=zeros(numel(rangesX)-1,numel(rangesY)-1);
    yPlotMean=zeros(numel(rangesX)-1,numel(rangesY)-1);
    xPlotMedian=zeros(numel(rangesX)-1,numel(rangesY)-1);
    yPlotMedian=zeros(numel(rangesX)-1,numel(rangesY)-1);
    zPlotMean=zeros(numel(rangesX)-1,numel(rangesY)-1);
    zPlotMedian=zeros(numel(rangesX)-1,numel(rangesY)-1);
    
    for ii=1:(numel(rangesX)-1)
        for jj=1:(numel(rangesY)-1)
            cellRanges{ii,jj}=find(x<rangesX(ii+1) & x>=rangesX(ii) ...
                & y<rangesY(jj+1) & y>=rangesY(jj));
            if numel(cellRanges{ii,jj})>0
                xPlotMean(ii,jj)=mean(x(cellRanges{ii,jj}));
                yPlotMean(ii,jj)=mean(y(cellRanges{ii,jj}));
                zPlotMean(ii,jj)=mean(z(cellRanges{ii,jj}));
                
                xPlotMedian(ii,jj)=median(x(cellRanges{ii,jj}));
                yPlotMedian(ii,jj)=median(y(cellRanges{ii,jj}));
                zPlotMedian(ii,jj)=median(z(cellRanges{ii,jj}));
            end
        end
    end
    
    rangesX(1)=rangesX(2);
    rangesY(1)=rangesY(2);
    rangesX(end)=rangesX(end-1);
    rangesY(end)=rangesY(end-1);
    for ii=1:(numel(rangesX)-1)
        for jj=1:(numel(rangesY)-1)
            if numel(cellRanges{ii,jj})==0
                
                xPlotMean(ii,jj)=mean(rangesX(ii:ii+1));
                yPlotMean(ii,jj)=mean(rangesY(jj:jj+1));
                
                xPlotMedian(ii,jj)=median(rangesX(ii:ii+1));
                yPlotMedian(ii,jj)=median(rangesY(jj:jj+1));
            end
        end
    end
    %l=plot(ax,xPlotMean,yPlotMean,xPlotMedian,yPlotMedian);
    l(1)=mesh(ax,xPlotMean,yPlotMean,zPlotMean);
    l(2)=mesh(ax,xPlotMedian,yPlotMedian,zPlotMedian);
    l(1).FaceColor='none';
    l(2).FaceColor='none';
end

function []=AddNumbersInRange(ax,overUnder,dim)
    hold(ax,'on')
    flag=true;
    l=findobj(ax,'type','line');
    if nargin<3
        dim=1;
    end
    fColor=[0.2 0.2 0.2];
    if isempty(l)
        l=findobj(ax,'type','scatter');
    end
    if dim==1
        x=[l.XData];
        y=[l.YData];
    else
        y=[l.XData];
        x=[l.YData];
    end
    over=sum(y>overUnder);
    under=numel(y(~isnan(y)))-over;
    
    boxAx=axis(ax);
    if dim==2
        boxAx=[boxAx(3:4),boxAx(1:2)];
    end
    minMax=[min(x(isfinite(x))),max(x(isfinite(x)))];
    boxAx(1:2)=minMax;
    if strcmp(ax.XScale,'linear')
        boxAx(1:2)=boxAx(1:2)+[-0.1 0.1]*(boxAx(2)-boxAx(1));
    else
        boxAx(1:2)=10.^(log10(boxAx(1:2))+([-0.1 0.1]*log10(boxAx(2)-boxAx(1))));
    end
    xDelta=boxAx(2)-0.02*(boxAx(2)-boxAx(1));
    yOver=(boxAx(4)-overUnder)*0.5+overUnder;
    yUnder=(overUnder-boxAx(3))*0.5+boxAx(3);
    if(yUnder)<boxAx(4) && under>0
        if dim==1
            text(ax,xDelta,yUnder,int2str(under),'HorizontalAlignment','right','color',fColor);
        else
            text(ax,yUnder,xDelta,int2str(under),'HorizontalAlignment','right','color',fColor);
        end
    else
        flag=false;
    end
    if(yUnder)>boxAx(3) && over>0
        if dim==1
            text(ax,xDelta,yOver,int2str(over),'HorizontalAlignment','right','color',fColor);
        else
            text(ax,yOver,xDelta,int2str(over),'HorizontalAlignment','right','color',fColor);
        end
    else
        flag=false;
    end
    if flag
        if dim==1
            plot(ax,minMax,[1 1]*overUnder,'--','color',fColor)
        else
            plot(ax,[1 1]*overUnder,minMax,'--','color',fColor)
        end
        
    end
    if dim==1
        axis(ax,boxAx)
    else
        axis(ax,boxAx([3:4,1:2]))
    end
    
    if dim==1
    else
    end
    
end



