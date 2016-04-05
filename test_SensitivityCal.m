global unstructglobal
sensSnax(:,find(sum(abs(sensSnax))==0))=[];
dCurr=[snaxel(:).d];
kk=1;
for ii=1:length(snaxel)
    snaxOrd(ii)=FindObjNum([],[snaxel(kk).snaxnext],[snaxel(:).index]);
    kk=snaxOrd(ii);
end
for ii=1:length(sensSnax(1,:)),
    snaxCopy=snaxel;
    e1=(1)./sensSnax(:,ii);
    e2=(-1)./sensSnax(:,ii);
    e1_sel=min(e1(e1>0));
    e2_sel=min(e2(e2>0));
    e_sel(ii)=min([e1_sel,e2_sel]);
    dChange{ii}=5*sensSnax(:,ii)/max(abs(sensSnax(:,ii)));
    dAct=dCurr'+dChange{ii};
    
    for jj=1:length(snaxel)
        snaxCopy(ii).d=dAct(ii);
    end
    [snakposition2]=PositionSnakes(snaxCopy,unstructglobal);
    %[snakposition2]=SnaxelNormal2(snaxCopy,snakposition2);
%     figh=CheckResultsLight(unstructglobal,snakposition,snaxel);
%     hold on
%     CheckResultsLight(unstructglobal,snakposition2,snaxCopy,figh);

coord1=vertcat(snakposition(:).coord);
coord2=vertcat(snakposition2(:).coord);
figure
plot(coord1(snaxOrd,1),coord1(snaxOrd,2),coord2(snaxOrd,1),coord2(snaxOrd,2))



end