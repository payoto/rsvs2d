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
    
    [desVarConstr,nMach]=ExtractVariables(varExtract,paramoptim);
    
    switch desVarConstr
        case 'MeanVolFrac'
            [loop(1)]=ConstantArea_Parabola(xMin,xMax,A,nPoints);
            [loop(2)]=ConstantArea_Wedge(xMin,xMax,A);
            resTag={'LinRes_ogive','LinRes_wedge'};
            
            parfor ii=1:length(loop)
                [obj(ii)]=OutputAndRunFlowSolve(loop(ii),rootFolder,tag{ii},paramoptim);
            end
            linTheoryOptim=([obj(:).cd]);
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
    
    [obj]=CutCellFlow_Handler(paramoptim,boundaryPath);
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

function [loop]=ConstantArea_Wedge(xMin,xMax,A)
    
    h=2*A/2/(xMax-xMin);
    
    points=[xMin,0;xMin+(xMax-xMin)/2,-h;xMax,0;xMin+(xMax-xMin)/2,h];
    
    loop.subdivision=points;
    loop.isccw=true;
    
end

function [loop]=ConstantArea_Klunker(xMin,xMax,A,M,nPoints)
    
    m = sqrt(M^2-1); %
    Pb = -4; % stagnation cp
    
    xl=(1-(m*Pb)/(12*A))/(1-(m*Pb)/(4*A));
    
    t=3/2*A*xl*(1-m*Pb/(12*A));
    
   xCoord=linspace(0,1,ceil(nPoints/2))';
   yCoord=t/2*xCoord/xl.*(2-xCoord/xl);
    
   points=[[xCoord,-yCoord];[xCoord(end:-1:2),yCoord(end:-1:2)]];
   
   ratio=xMax-xMin;
   points(:,1)=points(:,1)*ratio+xMin;
   points(:,2)=points(:,2)/ratio;
   
   loop.subdivision=points;
    loop.isccw=true;
end