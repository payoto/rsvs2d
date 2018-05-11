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
        ax22=subplot(1,2,2);
        hold on
        h3=figure('Name',['Convergence History eDV ',int2str(nIter(jj))]);
        for ii=1:4
            ax3(ii)=subplot(2,2,ii);
            hold on
        end
        h4=figure('Name',['Convergence History nonBasis design  ',int2str(nIter(jj))]);
        ax4=axes(h4);
        hold on
        
        c=get(ax2,'colororder');
        hold on
    
        clear l;
        kk=1;
        for iii=1:numel(lErrVec)
            ASOstruct=ASOstructAll(funcLogTest(lErrVec{iii}));
            ASOstruct=ASOstruct(find(~cellfun(@isempty,regexp({ASOstruct.location},...
                ['profile_',int2str(nIter(jj)),'$']))));
            if numel(lErrVec)>1
                fColor=c(mod(iii-1,size(c,1))+1,:);
            end

            for ii=1:numel(ASOstruct)
                nums=cellfun(@str2double,regexp(regexprep(ASOstruct(ii).location,...
                    '^.*profile_',''),'_','split'));
                l(kk)=plot(ax2,ASOstruct(ii).majorIt+ASOstruct(ii).DEIter,...
                    ASOstruct(ii).obj,'.-','color',fColor);
                plot(ax22,ASOstruct(ii).majorIt+ASOstruct(ii).DEIter,...
                    ASOstruct(ii).opt(ASOstruct(ii).majorIt),'.-','color',fColor);
                    %ones([1 numel(ASOstruct(ii).obj)])*nums(end),
                 l(kk).DisplayName=lName{iii};
                 kk=kk+1;
                 plot(ax4,ASOstruct(ii).nonBasisDesign,'.-','color',fColor)
                 for lll=1:4
                     plot(ax3(lll),ASOstruct(ii).eDV(:,lll),'.-','color',fColor)
                 end
            end
           
        end
        if exist('l','var')
            legend(ax2,l);
        end
    end
    
    
end

