%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2016
%
%          Supersonic Analytical Results
%
%             Alexandre Payot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [linTheoryOptim]=SupersonicOptimLinRes(paramoptim,rootFolder,xMin,xMax,A,nPoints)
    
    varExtract={'desVarConstr','nMach'};
    linTheoryOptim=0;
    [desVarConstr,nMach]=ExtractVariables(varExtract,paramoptim);
    for ii=1:length(desVarConstr)
        switch desVarConstr{ii}
            case 'MeanVolFrac'
                [loop(1)]=ConstantArea_Parabola(xMin,xMax,A,nPoints);
                [loop(2)]=ConstantArea_Wedge(xMin,xMax,A);
                [loop(3)]=ConstantArea_Klunker(xMin,xMax,A,nMach,nPoints);
                resTag={'LinRes_ogive','LinRes_wedge','LinRes_Klunker'};
                
                parfor jj=1:length(loop)
                    [obj(jj)]=OutputAndRunFlowSolve(loop(jj),rootFolder,resTag{jj},paramoptim);
                end
                linTheoryOptim=min([obj(:).cd]);
            case 'MinSumVolFrac'
                
                [loop(1)]=ConstantArea_Parabola(xMin,xMax,A,nPoints);
                [loop(2)]=ConstantArea_Wedge(xMin,xMax,A);
                [loop(3)]=ConstantArea_Klunker(xMin,xMax,A,nMach,nPoints);
                resTag={'LinRes_ogive','LinRes_wedge','LinRes_Klunker'};
                
                parfor jj=1:length(loop)
                    [obj(jj)]=OutputAndRunFlowSolve(loop(jj),rootFolder,resTag{jj},paramoptim);
                end
                linTheoryOptim=min([obj(:).cd]);
                
        end
    end
    
    %[loop]=ConstantArea_Klunker(xMin,xMax,A,nMach,nPoints);
    
    
end

function [obj]=OutputAndRunFlowSolve(loop,rootFolder,tag,paramoptim)
    
    boundaryFolder=[rootFolder,filesep,tag];
    mkdir(boundaryFolder);
    boundaryPath=[boundaryFolder,filesep,'boundary_',tag,'.dat'];
    fid=fopen(boundaryPath,'w');
    loop.subdivision(end+1,:)=loop.subdivision(1,:);
    BoundaryOutput(loop,fid);
    
    [obj]=CutCellFlow_Handler(paramoptim,boundaryFolder);
end

%% Area Constraint - Profile Generation

function [loop]=ConstantArea_Parabola(xMin,xMax,A,nPoints)
    
    f=@(x) [x.^2, x, ones(size(x))];
    F=@(x) [(x.^3)/3, (x.^2)/2, x ];
    
    R=[F(xMax)-F(xMin);f(xMax);f(xMin)];
    targ=[A/2;0;0];
    
    coeff=R\targ;
    
    xCoord=linspace(xMin,xMax,ceil(nPoints/2))';
    yCoord=f(xCoord)*coeff;
    
    points=[[xCoord,-yCoord];[xCoord(end-1:-1:2),yCoord(end-1:-1:2)]];
    
    loop.subdivision=points;
    loop.isccw=true;
    
end

function [loop]=ConstantArea_Wedge(xMin,xMax,A,M)
    
    h=2*A/2/(xMax-xMin);
    Dx=xMax-xMin;
    tanDel=h/Dx;
    
    points=[xMin,0;xMin+(xMax-xMin)/2,-h;xMax,0;xMin+(xMax-xMin)/2,h];
    
    loop.subdivision=points;
    loop.isccw=true;
    
end


function [loop]=ConstantArea_Busemann(xMin,xMax,A)
    
    h=A/2/(xMax-xMin);
    
    points=[xMin,0;xMax,0;xMin+(xMax-xMin)/2,h];
    
    loop.subdivision=points;
    loop.isccw=true;
    
end

function [tanDel]=CalcTanDel(M,B)
    
    tanDel=2/tan(B)*(M^2*(sin(B)^2)-1)/(M^2*(1.4+cos(2*B))+2);
    
end

function [loop]=ConstantArea_Klunker(xMin,xMax,A,M,nPoints)
    ratio=xMax-xMin;
    
    Acalc=A;
    
    m = sqrt(M^2-1); %
    Pb = 0.7/-(M^2-1); % base cp
    
    xl=(1-(m*Pb)/(12*Acalc))/(1-(m*Pb)/(4*Acalc));
    
    t=3/2*Acalc*xl*(1-m*Pb/(12*Acalc));
    
   xCoord=linspace(0,1,ceil(nPoints/2))';
   yCoord=t/2*xCoord/xl.*(2-xCoord/xl);
    
   points=[[xCoord,-yCoord];[xCoord(end:-1:2),yCoord(end:-1:2)]];
   
   ratio=xMax-xMin;
   points(:,1)=points(:,1)*ratio+xMin;
   points(:,2)=points(:,2)/ratio;
   
   [A2]=CalculatePolyArea(points);
   loop.subdivision=points;
    loop.isccw=true;
end