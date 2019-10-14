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


function [linTheoryOptim,obj]=SupersonicOptimLinRes(paramoptim,rootFolder,xMin,xMax,A,nPoints)
    
    varExtract={'desVarConstr','nMach'};
    linTheoryOptim=0;
    [desVarConstr,nMach]=ExtractVariables(varExtract,paramoptim);
    for ii=1:length(desVarConstr)
        switch desVarConstr{ii}
            case 'MeanVolFrac'
                [loop{1}]=ConstantArea_Parabola(xMin,xMax,A,nPoints);
                [loop{2}]=ConstantArea_Wedge(xMin,xMax,A);
                [loop{3}]=ConstantArea_Klunker(xMin,xMax,A,nMach,nPoints);
                [loop{4}]=ConstantArea_Busemann(xMin,xMax,A,nMach);
                resTag={'LinRes_ogive','LinRes_wedge','LinRes_Klunker','LinRes_Busemann'};
                
                parfor jj=1:length(loop)
                    [obj(jj)]=OutputAndRunFlowSolve(loop{jj},rootFolder,resTag{jj},paramoptim);
                end
                linTheoryOptim=min([obj(:).cd]);
            case 'ValVolFrac'
                [loop{1}]=ConstantArea_Parabola(xMin,xMax,A,nPoints);
                [loop{2}]=ConstantArea_Wedge(xMin,xMax,A);
                [loop{3}]=ConstantArea_Klunker(xMin,xMax,A,nMach,nPoints);
                [loop{4}]=ConstantArea_Busemann(xMin,xMax,A,nMach);
                resTag={'LinRes_ogive','LinRes_wedge','LinRes_Klunker','LinRes_Busemann'};
                
                parfor jj=1:length(loop)
                    [obj(jj)]=OutputAndRunFlowSolve(loop{jj},rootFolder,resTag{jj},paramoptim);
                end
                linTheoryOptim=min([obj(:).cd]);
            case 'MinSumVolFrac'
                
                [loop{1}]=ConstantArea_Parabola(xMin,xMax,A,nPoints);
                [loop{2}]=ConstantArea_Wedge(xMin,xMax,A);
                [loop{3}]=ConstantArea_Klunker(xMin,xMax,A,nMach,nPoints);
                [loop{4}]=ConstantArea_Busemann(xMin,xMax,A,nMach);
                resTag={'LinRes_ogive','LinRes_wedge','LinRes_Klunker','LinRes_Busemann'};
                
                parfor jj=1:length(loop)
                    [obj(jj)]=OutputAndRunFlowSolve(loop{jj},rootFolder,resTag{jj},paramoptim);
                end
                
                linTheoryOptim=min([obj(:).cd]);
                MinValVolFrac
            case 'MinValVolFrac'
                
                [loop{1}]=ConstantArea_Parabola(xMin,xMax,A,nPoints);
                [loop{2}]=ConstantArea_Wedge(xMin,xMax,A);
                [loop{3}]=ConstantArea_Klunker(xMin,xMax,A,nMach,nPoints);
                [loop{4}]=ConstantArea_Busemann(xMin,xMax,A,nMach);
                resTag={'LinRes_ogive','LinRes_wedge','LinRes_Klunker','LinRes_Busemann'};
                
                parfor jj=1:length(loop)
                    [obj(jj)]=OutputAndRunFlowSolve(loop{jj},rootFolder,resTag{jj},paramoptim);
                end
                
                linTheoryOptim=min([obj(:).cd]);
        end
    end
    
    %[loop]=ConstantArea_Klunker(xMin,xMax,A,nMach,nPoints);
    
    
end
