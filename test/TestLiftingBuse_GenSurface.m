function []=TestLiftingBuse_GenSurface(teststruct, testCellParam)
    kk=0;
    listNames=cell(0);
    
    testCellFunc = @(testCellParam) (numel(testCellParam{1,2})*numel(testCellParam{2,2}));
    
    for ii=1:(numel(teststruct)/testCellFunc(testCellParam))
        indStart=(ii-1)*testCellFunc(testCellParam)+1;
        indEnd=(ii)*testCellFunc(testCellParam);
        inds = indStart:indEnd;
        inds = reshape(inds, [numel(testCellParam{1,2}),numel(testCellParam{2,2})]);
        seriesName=[regexprep(regexprep(teststruct(indStart).geometry,...
            '^.*boundary_',''),'\.dat',''),...
            '_surface'];
        
        jj=0;
        flagMatch=false;
        while jj<numel(listNames) && ~flagMatch
            jj=jj+1;
            flagMatch=strcmp(listNames{jj},seriesName);
        end
        
        if ~flagMatch
            listNames{end+1}=seriesName;
            h(jj+1)=figure('Name',seriesName);
        else
            figure(jj);
        end
        fieldsRes=fieldnames(teststruct(indStart).res);
        for kk=1:numel(fieldsRes)
            subplot(3,3,kk)
            x=zeros(size(inds));
            y=zeros(size(inds));
            z=zeros(size(inds));
            for ll=1:numel(inds)
                x(ll)=[teststruct(inds(ll)).nAlpha];
                y(ll)=[teststruct(inds(ll)).area];
                z(ll)=[teststruct(inds(ll)).res.(fieldsRes{kk})];
            end
            surf(x,y,z,'DisplayName',[teststruct(indStart).mesher,'_',...
                num2str(teststruct(indStart).area),'_',int2str(kk)])
            title(fieldsRes{kk})
            zlabel(fieldsRes{kk})
            xlabel('incidence')
            ylabel('area')
            hold on
        end
        if numel(fieldsRes)< 9
            subplot(3,3,numel(fieldsRes)+1)
            cl=zeros(size(inds));
            cd=zeros(size(inds));
            for ll=1:numel(inds)
                cl(ll)=[teststruct(inds(ll)).res.cl];
                cd(ll)=[teststruct(inds(ll)).res.cd];
            end
            surf(y,cd,cl,'DisplayName',[teststruct(indStart).mesher,'_',...
                num2str(teststruct(indStart).area),'_',int2str(kk)])
            title('drag polar')
            zlabel('cl')
            xlabel('area')
            ylabel('Cd')
            hold on
        end
    end
end