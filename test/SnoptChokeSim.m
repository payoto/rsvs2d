
function [saveX] = SnoptChokeSim(startPt, isConstrAct, isNewPlot, is3D, ax, plotFormat)
    global SnoptChokeSim_HASRUNBEFORE;
    addpath ../../ASO_LK/SRC/matlab-snopt/
    addpath ../../ASO_LK/SRC/
    if isempty(SnoptChokeSim_HASRUNBEFORE)
        SnoptChokeSim_HASRUNBEFORE=false;
    end
    discontinuityMulti = @(x) ((tanh(100*(1-sum(x.^2,1)))+1)/2);
%     discontinuityMulti = @(x) (sum(x.^2,1)<1);
    userfun = @(x) ChokeBoundModelFunctionSNOPT(x);
    
    nDim = 2;

    dim = ones(nDim,1);
    
    if nargin == 0;
        startPt = 0*dim;
    elseif numel(startPt)==1;
        startPt = (startPt) * dim;
    end
    if nargin<6
        plotFormat = {'linewidth',2,'marker','o'};
    end
    if nargin<4
        is3D = 0;
    end
    if nargin<2 % activate circular constraint which runs along the discontinuity
        isConstrAct = 0;
    end
    if nargin<3
        if SnoptChokeSim_HASRUNBEFORE
            isNewPlot=0;
            
        else
            isNewPlot=1;
        end
    end
    
    if isConstrAct==1
        plotFormat(end+1:end+2) = {'linestyle','--'};
    end

    saveX = zeros([nDim+1,0]);
    saveLSX = zeros([nDim+1,0]);
   
    plotPoints = @(points,form)plot(points([1:end],1),points([1:end],2),form{:});
   
    startPt = reshape(startPt,size(dim));
    x = startPt;
    xlow = -10 * dim;
    xupp = 10 * dim;
    xmul = 1 * dim;
    xstate = 1 * dim;

    Flow= [-inf];
    Fupp= [inf];
    Fmul=[1];
    Fstate=[1];

    A=[];
    iAfun=[];
    jAvar=[];
    iGfun=[]; 
    jGvar=[];

    % [x,F,inform,xmul,Fmul] = snopt( x, xlow, xupp, xmul, xstate,...
    % 					 Flow, Fupp, Fmul, Fstate,...
    % 					 userfun,...
    %                      A, iAfun, jAvar, iGfun, jGvar...
    %                  );

    optFun = userfun;
    x0 = x;
    A = [];
    b= [] ;
    Aeq = [];
    beq = [];
    lb = xlow;
    ub = xupp;

    nlconFun = @nlcon;

    [xEnd,fval,exitflag,lambda,nFunCall] = snfmincon(optFun, x0, A, b, Aeq,beq,...
        lb,ub,nlconFun...
        );
    plotPoints = @(ax,points,form)plot(ax,points([1:end],1),points([1:end],2),form{:});
    plotPoints3 = @(ax,points,form)plot3(ax,points([1:end],1),...
        points([1:end],2),points([1:end],3),form{:});
    if isNewPlot
        if ~exist('ax','var')
            h = figure;
            hold on;
        end
        
        [x,y] = meshgrid(unique(linspace(-2,2,100)));
        z = reshape(ChokeBoundModelFunctionNoGrad([x(:),y(:)]),size(x));
        if ~is3D
            contour(x,y,z,30)
        else
            surfc(x,y,z)
        end
    end
    if ~exist('ax','var')
        ax = gca;
    else
        ax.NextPlot = 'add';
        
    end
    if ~is3D
%         l = plotPoints(ax,saveLSX', {'k--'});
        l = plotPoints(ax,saveX', plotFormat);
    else
%         l = plotPoints3(ax,saveLSX', {'k--'});
        l = plotPoints3(ax,saveX', plotFormat);
    end
    if ~is3D
        plot(ax,xEnd(1),xEnd(2),'Color', l.Color,'markersize',10,'marker','x', 'linewidth',2);
    else
        box = axis;
        plot3(ax,xEnd(1),xEnd(2),box(5),'Color', l.Color,'markersize',10,'marker','x', 'linewidth',2);

    end
    constrStr = 'noconstr';
    if isnan(isConstrAct)
        constrStr = 'fail';
    elseif (isConstrAct)
        constrStr = 'constr';
    end
        
    l.DisplayName = [...
        sprintf('$(%.2f,\\ %.2f)\\rightarrow(%.3f,\\ %.3f)$'')'...
        ,startPt(1),startPt(2),xEnd(1),xEnd(2)),...
        '; ',...
        constrStr, '; $f^*=', num2str(fval),'$'];
    ls = findobj(ax,'type','line');
    leg = legend(ls(~cellfun(@isempty,{ls.DisplayName})));
    leg.Interpreter = 'latex';
    if ~is3D
        axis([-2 2 -2 2])
    else
        box = axis;
        box(1:4) = [-2 2 -2 2];
        axis(box)
    end

    function [out, G]=ChokeBoundModelFunctionSNOPT(x)
        
        f = @(x) 0.5*sum(x,1) + sum(x.^2,1) + discontinuityMulti(x).*(10- 3*sum(x.^2,1));
        
        out = f(x);
        
        G = zeros(size(x));
        eps = 1e-8;
        deltaX = zeros(size(x));
        for ii = 1:numel(size(x))
            deltaX(ii) = eps;
            G(ii) = (f(x+deltaX) - f(x-deltaX))/(2*eps);
            deltaX(ii)=0;
        end
        if isnan(isConstrAct) && (sum(x.^2)<1)
            out=nan;
            G=G*nan;
        end
        
        
        x = [x;out];
        if numel(saveX)==0 || out<=saveX(end,end)
            saveX(:,end+1)=x;
        else
            saveLSX(:,end+1) = x;
        end
    end
    
    function [c, ceq, gc, gceq] = nlcon(x)
        if ~isnan(isConstrAct)
            c = -sum(x.^2)+isConstrAct*2-1;
        elseif -sum(x.^2)+1 > 0;
            c = nan;
            x=nan;
        else
            c = -sum(x.^2)-1;
        end
        ceq = [];
        gc = -2.*x;
        gceq = [];
    end
    function [out]=ChokeBoundModelFunctionNoGrad(x)

            f = @(x) 0.5*sum(x,2) + sum(x.^2,2) + discontinuityMulti(x')'.*(10- 3*sum(x.^2,2));
            
            out = f(x);
            if isnan(isConstrAct)
                out((sum(x.^2,2)<1))=nan;
            end

    end
    SnoptChokeSim_HASRUNBEFORE = true;
end







