
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
    nVar=11;
    nNewPoints=2;
    coeffVec=sym('a',[nUniqVar,1],'real');
    coeff2=sym(zeros([2*nVar,1]));
    if mod(nVar,2)==0
        newStencil=[coeffVec(end:-1:1);coeffVec];
        coeff2(1:2:end)=newStencil;
    else
        newStencil=[coeffVec(end:-1:2);coeffVec];
        coeff2(1:2:end)=newStencil;
    end
    %coeff2=coeffVec
    numPoints=2*n;
    numNewPoints=nNewPoints*2*n;
    subMask=zeros(numNewPoints,numPoints);
    subMask=sym(subMask);
    for ii=0:n-1
        for kk=1:2
            iStart=2*ii+(kk-1);
            jStart=2*ii*nNewPoints+kk-1-nVar+3-1+mod(nVar,2);

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

% sol=solve([coeffsolMat(:,1)==areaMat(:,1); sum(subMask,2)==1; subMask(3,1)~=1  ],'ReturnConditions',true,'Real',true)
sol=solve([coeffsolMat(:,1)==areaMat(:,1); sum(subMask,2)==1  ],'ReturnConditions',true,'Real',true) % not tested precendent line has been used

%% Plot all solutions

solVec=sym('a',[nUniqVar,1],'real');
subMask2=subMask(1:2:end,1:2:end);
[~,~,solNum]=unique(sol.conditions);
h = figure;
%points=[0 0;1 0;2 0;2 2;1 2; 1 1;0.5 1;0 1];
points=[0 0;1 0; 1 1;0 1];

axisLims = [min(points)-1,max(points)+1];
axisLims = axisLims([1,3,2,4]);

for solI=2:numel(sol.conditions);
    solI
    ax(solI)=subplot(1,numel(sol.conditions),solI);
    % solI=4;
    %isAlways(subs(sol.conditions, 'x', -1))
    paramsAct=1:numel(sol.parameters);
    if  nUniqVar==6 && nVar==12
        switch solNum(solI)
            case 1 % x < 0
                parameterTest = -1./[0.01 0.1 0.5 1 2 4 8 16 32 64 128 1000];
            case 2 % 1 <= x
                parameterTest = 1 + 1./[inf 100 10 5 2 1 0.5 0.25];
            case 3 % 1/2 <= x
                parameterTest = 0.5 + 1./[inf 100 10 5 2 1 0.5 0.25];
        end
    else
        parameterTest=(lhsdesign(20,4)'-0.5)*2;
    end
    parameterTest = 0;
    newStencilAct=newStencil;
    for ii=1:nUniqVar

        solVec(ii)=eval(['sol.a',int2str(ii),'(',int2str(solI),')']);
        subMask2=subs(subMask2,['a',int2str(ii)],solVec(ii));
        newStencilAct=subs(newStencilAct,['a',int2str(ii)],solVec(ii));
    end

    solVec
    for x=parameterTest
%         valStenc=subs(newStencilAct,sol.parameters(paramsAct),x');
        valStenc = newStencilAct;
        strct.varStencil=double(vpa(valStenc));
        strct.nNew=nNewPoints;
        [newPoints]=SubSurfVarStencil_NoCorn_STA(points,2,strct);

%         figure(hMain);
        plot(ax(solI),newPoints(:,1),newPoints(:,2), 'DisplayName',num2str(x));
        hold on
%         test_subdiv_varstenc(points, strct)

    end
    axis(axisLims)
end

TightenAxes(ax,'auto');
for ii = 1:numel(ax)
    ax(ii).XLabel.String = (['Solution ', int2str(ii)]);
end
    


%% Choose solution 5 with ternary

solI=5;

solVec=sym('a',[nUniqVar,1],'real');
subMask2=subMask(1:2:end,1:2:end);
[~,~,solNum]=unique(sol.conditions);
h = figure;
axSol5 = subplot(121);
hold on;
axSol5_2 = subplot(122);
hold on;
points=[0 0;1 0;2 0;2 2;1 2; 1 1;2/3 1 ; 1/3 1;0 1];
% points=[0 0;1 0; 1 1;0 1];



paramsAct=1;
newStencilAct=newStencil;
for ii=1:nUniqVar

    solVec(ii)=eval(['sol.a',int2str(ii),'(',int2str(solI),')']);
    subMask2=subs(subMask2,['a',int2str(ii)],solVec(ii));
    newStencilAct=subs(newStencilAct,['a',int2str(ii)],solVec(ii));
end
paramPrecusror = [41 43 42];%linspace(38,   50,10);
parameterTest = -1./paramPrecusror;
kk=1;
clear lineSol5
clear lineSol52
clear newfig
figRep = figure;
A = [];
Asubdiv = [];
for x=parameterTest
    valStenc=subs(newStencilAct,sol.parameters(paramsAct),x);
    strct.varStencil=double(vpa(valStenc));
    
    strct.nNew=nNewPoints;
    A(kk)= CalculatePolyArea(points);
    [newPoints]=SubSurfVarStencil_NoCorn_STA(points,4,strct);
    Asubdiv(kk) = CalculatePolyArea(newPoints);


%         figure(hMain);
    lineSol5(kk)=plot(axSol5,newPoints(:,1),newPoints(:,2), 'DisplayName',num2str(x));
    hold on
	[figRep, newfig(kk), linesCurv]=test_subdiv_varstenc(points, strct,7,figRep,true);
    [linesCurv.DisplayName] = deal(num2str(paramPrecusror(kk)));
    lineSol52(kk) = linesCurv(end);
    
    kk = kk+1;
end
legend(lineSol5)

legend(lineSol52)
% Eplore the properties of the stencil itself
kk=1;
delete(findobj(axSol5_2,'type','line'))
varStencil = zeros(size(strct.varStencil,1),numel(parameterTest));
for x=parameterTest
    valStenc=subs(newStencilAct,sol.parameters(paramsAct),x);
    varStencil(:,kk)=double(vpa(valStenc));
    kk=kk+1;
end

plot(axSol5_2,parameterTest,varStencil(size(varStencil,1)/2+1:end,:)')

A
Asubdiv