
clc
close all
clear all

%%
points=[0 0;1 0;2 0;2 2;1 2; 1 1;0.5 1;0 1]
pointsVec=points';
pointsVec=pointsVec(:);
plot(points(:,1),points(:,2));
n=length(points(:,1));
centreMat=eye(2*n);
centreMat=(centreMat+centreMat(:,[end-1:end,1:end-2]))*0.5;

[rotDif]=[0 -1 0 1; 1 0 -1 0];
normMat=zeros(2*n);
for ii=1:n-1
    normMat((2*(ii-1)+1):(2*(ii-1)+2),(2*(ii-1)+1):(2*(ii-1)+4))=rotDif;
end
ii=n;
normMat((2*(ii-1)+1):(2*(ii-1)+2),(2*(ii-1)+1):(2*(ii-1)+2))=rotDif(:,1:2);
normMat((2*(ii-1)+1):(2*(ii-1)+2),1:2)=rotDif(:,3:4);
A=0.5*(normMat*pointsVec)'*(centreMat*pointsVec);
A=0.5*pointsVec'*(normMat'*centreMat)*pointsVec;

areaMat=(normMat'*centreMat);

%% even
isEven=false;
if isEven
    nVar=5;
    nUniqVar=nVar;
    nNewPoints=2;
    coeffVec=sym('a',[nVar,1],'real');
    coeff2=sym(zeros([4*nVar,1]));
    
    newStencil=[coeffVec(nVar:-1:1);coeffVec];
    coeff2(1:2:end)=newStencil;
    %coeff2=coeffVec
    numPoints=2*n;
        numNewPoints=nNewPoints*2*n;
    subMask=zeros(numNewPoints,numPoints);
    subMask=sym(subMask);
    for ii=0:n-1
        for kk=1:2
            iStart=2*ii+(kk-1);
            jStart=2*ii*nNewPoints+kk-1-2*nVar+nNewPoints;

            nJ=4*nVar;
            indX=zeros(1,1);
            indY=zeros(1,nJ);
            for iLoop=1:1
                indX(iLoop)=mod(iStart+(iLoop-1),numPoints)+1;
            end
            for jLoop=1:nJ
                indY(jLoop)=mod(jStart+(jLoop-1),numNewPoints)+1;
            end

            subMask(indY,indX)=(coeff2)+(subMask(indY,indX));
        end
    end
else
    nUniqVar=6;
    nVar=12;
    nNewPoints=3;
    coeffVec=sym('a',[nUniqVar,1],'real');
    coeff2=sym(zeros([2*nVar,1]));
    newStencil=[coeffVec(end:-1:1);coeffVec]
    coeff2(1:2:end)=newStencil;
    %coeff2=coeffVec
    numPoints=2*n;
    numNewPoints=nNewPoints*2*n;
    subMask=zeros(numNewPoints,numPoints);
    subMask=sym(subMask);
    for ii=0:n-1
        for kk=1:2
            iStart=2*ii+(kk-1);
            jStart=2*ii*nNewPoints+kk-1-nVar+3-1;

            nJ=length(coeff2);
            indX=zeros(1,1);
            indY=zeros(1,nJ);
            for iLoop=1:1
                indX(iLoop)=mod(iStart+(iLoop-1),numPoints)+1;
            end
            for jLoop=1:nJ
                indY(jLoop)=mod(jStart+(jLoop-1),numNewPoints)+1;
            end

            subMask(indY,indX)=(coeff2)+(subMask(indY,indX));
        end
    end
    
end
subMask(1:2:end,1:2:end)
%%

n2=nNewPoints*length(points(:,1));
centreMatLarge=eye(2*n2);
centreMatLarge=(centreMatLarge+centreMatLarge(:,[end-1:end,1:end-2]))*0.5;

[rotDif]=[0 -1 0 1; 1 0 -1 0];
normMatLarge=zeros(2*n2);
for ii=1:n2-1
    normMatLarge((2*(ii-1)+1):(2*(ii-1)+2),(2*(ii-1)+1):(2*(ii-1)+4))=rotDif;
end
ii=n2;
normMatLarge((2*(ii-1)+1):(2*(ii-1)+2),(2*(ii-1)+1):(2*(ii-1)+2))=rotDif(:,1:2);
normMatLarge((2*(ii-1)+1):(2*(ii-1)+2),1:2)=rotDif(:,3:4);


areaMatLarge=(normMatLarge'*centreMatLarge);

areaSym=sym(areaMatLarge);

coeffsolMat=subMask'*areaMatLarge*subMask;

sol=solve([coeffsolMat(:,1)==areaMat(:,1); sum(subMask,2)==1; subMask(3,1)~=1  ],'ReturnConditions',true,'Real',true)

%% Exploit sol

solVec=sym('a',[nUniqVar,1],'real');
subMask2=subMask(1:2:end,1:2:end);
%for solI=1:5;
solI=5;
paramsAct=1;
 newStencilAct=newStencil;
for ii=1:nUniqVar
    
    solVec(ii)=eval(['sol.a',int2str(ii),'(',int2str(solI),')']);
    subMask2=subs(subMask2,['a',int2str(ii)],solVec(ii));
    newStencilAct=subs(newStencilAct,['a',int2str(ii)],solVec(ii));
end

%hMain=figure;
hMain=1;
solVec
figKeep = figure;
for ii=linspace(36,40,10)
    %points=[0 0;1 0;2 0;2 2;1 2; 1 1;0.5 1;0 1];
    points=[0 0;1 0; 1 1;0 1];
    valStenc=subs(newStencilAct,sol.parameters(paramsAct),[-1/ii])
    strct.varStencil=double(vpa(valStenc));
    strct.nNew=nNewPoints;
    [newPoints]=SubSurfVarStencil_NoCorn_STA(points,8,strct);
    
    figure(hMain);
    plot(newPoints(:,1),newPoints(:,2))
    hold on
    test_subdiv_varstenc(points, strct, 7,figKeep );

end

figure(hMain);
axis equal
plot(points(:,1),points(:,2));

%%
vecX=linspace(-1,0.0001,50);
vecY=linspace(-1,0.0001,50);
points=oldPoints;
[X,Y]=meshgrid(vecX,vecY);
%for mm=0.05:0.1:0.5
surfDat=zeros([length(vecX)*[1 1],nUniqVar]);
kk=1;
for ii=vecX
    ll=1;
    for jj=vecY
        surfDat(kk,ll,:)=double(subs(solVec,sol.parameters(paramsAct),[ii]));
        ll=ll+1;
    end
    kk=kk+1;
end

figure
for ii=1:nUniqVar
    subplot(2,4,ii)
    surf(X,Y,surfDat(:,:,ii));
    
end
%%
figure
contList={[0.5:0.05:1],[0.3:0.05:1],[0:0.05:0.5],[0:0.025:0.25],[-0.25:0.05:0],[-0.25:0.1:0.25]};
for ii=1:nUniqVar
    axh(ii)=subplot(2,3,ii);
    [c{ii},hMain{ii}]=contourf(X,Y,real(surfDat(:,:,ii)));
    set(hMain{ii},'levellist',contList{ii})
    clabel(c{ii},hMain{ii})
    
end
linkaxes(axh,'xy')
%end
%valx=subs(solVec,sol.parameters,[0.75 0.25 0.125])
