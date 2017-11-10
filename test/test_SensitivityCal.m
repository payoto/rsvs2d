global unstructglobal
sensSnax(:,find(sum(abs(sensSnax))==0))=[];
dCurr=[snaxel(:).d];
kk=1;
for ii=1:length(snaxel)
    snaxOrd(ii)=FindObjNum([],[snaxel(kk).snaxnext],[snaxel(:).index]);
    kk=snaxOrd(ii);
end
snaxOrd(end+1)=snaxOrd(1);
coord1=vertcat(snakposition(:).coord);
dir=vertcat(snakposition(:).vector);
testPos{2}=[];
for jj=1:2
    testPos{jj}=zeros(length(snaxel));
    for ii=1:length(snaxel)
        testPos{jj}(:,ii)=(coord1(ii,jj)>coord1(:,jj));
    end
end
for ii=1:length(sensSnax(1,:)),
    snaxCopy=snaxel;
    e1=(1)./sensSnax(:,ii);
    e2=(-1)./sensSnax(:,ii);
    e1_sel=min(e1(e1>0));
    e2_sel=min(e2(e2>0));
    e_sel(ii)=min([e1_sel,e2_sel]);
    dChange{ii}=sensSnax(:,ii)/max(abs(sensSnax(:,ii)));
    %dChange{ii}=-sensSnax(:,ii)/100;
    dAct=dCurr'+dChange{ii};
    
    for jj=1:length(snaxel)
        snaxCopy(jj).d=dAct(jj);
    end
    [snakposition2]=PositionSnakes(snaxCopy,unstructglobal);
    %[snakposition2]=SnaxelNormal2(snaxCopy,snakposition2);
    %     figh=CheckResultsLight(unstructglobal,snakposition,snaxel);
    %     hold on
    %     CheckResultsLight(unstructglobal,snakposition2,snaxCopy,figh);
    
    
    coord2=vertcat(snakposition2(:).coord);
    for jj=1:2
        testPos2{jj}=zeros(length(snaxel));
        for kk=1:length(snaxel)
            testPos2{jj}(:,kk)=(coord2(kk,jj)>coord2(:,jj));
        end
    end
    figure
    plot(coord1(snaxOrd,1),coord1(snaxOrd,2),'o-',coord2(snaxOrd,1),coord2(snaxOrd,2),'o-')
    hold on
    for jj=1:length(coord1(:,1))
        plot([coord1(jj,1),coord2(jj,1)],[coord1(jj,2),coord2(jj,2)],'k--')
    end
    title(['mode ',int2str(ii)])
    axis equal
    
end

