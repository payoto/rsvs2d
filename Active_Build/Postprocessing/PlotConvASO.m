function []=PlotConvASO(ASOstruct,nIter)
    
    
    
    ASOstructAll=ASOstruct;
    splitCase='RunName';
    switch splitCase
        case 'errorVecMode'
            errVecModes={ASOstruct.errorVecMode};
            lErrVec=unique(errVecModes);
            funcLogTest=@(in) ~cellfun(@isempty,regexp(errVecModes,in));
        case 'subDivLevel'
            for ii=1:numel(ASOstruct)
                errVecModes{ii}=int2str(unique(ASOstruct(ii).refLvl));
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
    
    fColor=[ 0 0 0];
    lName=regexprep(lErrVec,'_',' ');
    for jj=1:numel(nIter)
        h2=figure('Name',['Convergence History profile ',int2str(nIter(jj))]);
        ax2=subplot(1,2,1);
        hold on
        ax22(1)=subplot(2,2,2);
        hold on
        ax22(2)=subplot(2,2,4);
        hold on
        h3=figure('Name',['Convergence History eDV ',int2str(nIter(jj))]);
        for ii=1:4
            ax3(ii)=subplot(2,2,ii);
            hold on
        end
        h4=figure('Name',['Convergence History nonBasis design  ',int2str(nIter(jj))]);
        ax4=axes(h4);
        hold on
        h5=figure('Name',['Final Profiles  ',int2str(nIter(jj))]);
        ax5=axes(h5);
        hold on
        
        c=get(ax2,'colororder');
        markerOrd='.+d';
        hold on
    
        clear l;
        kk=1;
        for iii=1:numel(lErrVec)
            ASOstruct=ASOstructAll(funcLogTest(lErrVec{iii}));
            ASOstruct=ASOstruct(find(~cellfun(@isempty,regexp({ASOstruct.location},...
                ['profile_',int2str(nIter(jj)),'$']))));
            if numel(lErrVec)>0
%                 fColor=c(mod(iii-1,size(c,1))+1,:);
%                 fMarker=markerOrd(mod(floor((iii-1)/size(c,1)),numel(markerOrd))+1);
                fColor=c(mod(floor((iii-1)/numel(markerOrd)),size(c,1))+1,:);
                fMarker=markerOrd(mod((iii-1),numel(markerOrd))+1);
            end

            for ii=1:numel(ASOstruct)
                nums=cellfun(@str2double,regexp(regexprep(ASOstruct(ii).location,...
                    '^.*profile_',''),'_','split'));
                l(kk)=plot(ax2,ASOstruct(ii).majorIt+ASOstruct(ii).DEIter,...
                    ASOstruct(ii).obj,[fMarker,'-'],'color',fColor);
                for ll=1:numel(ASOstruct(ii).loops{end})
                    plot(ax5,ASOstruct(ii).loops{end}{ll}(:,1),...
                        ASOstruct(ii).loops{end}{ll}(:,2),[fMarker,'-'],'color',fColor);
                end
%                 plot(ax22(1),ASOstruct(ii).majorIt+ASOstruct(ii).DEIter,...
%                     ASOstruct(ii).opt(ASOstruct(ii).majorIt),[fMarker,'-'],'color',fColor);
                
                field={'adjointMax','adjoint'};
                for jjj=1:numel(field)
                    try
                [ASOstruct(ii).residual(cellfun(@isempty,...
                    {ASOstruct(ii).residual.(field{jjj})})).(field{jjj})]=deal(nan);
                res=[ASOstruct(ii).residual.(field{jjj})];
                plot(ax22(mod(jjj,numel(ax22))+1),ASOstruct(ii).majorIt+ASOstruct(ii).DEIter,...
                    res(ASOstruct(ii).majorIt),[fMarker,'-'],'color',fColor);
                    catch
                        warning(['field ', field{jj},' failed'])
                    end
                end
                    %ones([1 numel(ASOstruct(ii).obj)])*nums(end),
                 l(kk).DisplayName=lName{iii};
                 kk=kk+1;
                 plot(ax4,ASOstruct(ii).nonBasisDesign,[fMarker,'-'],'color',fColor)
                 for lll=1:4
                     plot(ax3(lll),ASOstruct(ii).eDV(:,lll),[fMarker,'-'],'color',fColor)
                 end
            end
           
        end
        if exist('l','var')
            legend(ax2,l);
        end
    end
    
    
end

