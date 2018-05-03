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
        ax2=axes(h2);
        
        c=get(ax2,'colororder');
        hold on
    
        clear l;
        kk
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
                l(iii)=plot(ax2,ASOstruct(ii).majorIt+ASOstruct(ii).DEIter,...
                    ASOstruct(ii).obj,'.-','color',fColor);
                    %ones([1 numel(ASOstruct(ii).obj)])*nums(end),
                 l(iii).DisplayName=lName{iii};
            end
           
        end
        if exist('l','var')
            legend(ax2,l);
        end
    end
    
    
end