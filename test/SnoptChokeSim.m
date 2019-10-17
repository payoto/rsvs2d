
function [saveX] = SnoptChokeSim(startPt, isConstrAct, isNewPlot, plotFormat)
    global SnoptChokeSim_HASRUNBEFORE;
    addpath ../../ASO_LK/SRC/matlab-snopt/
    addpath ../../ASO_LK/SRC/
    if isempty(SnoptChokeSim_HASRUNBEFORE)
        SnoptChokeSim_HASRUNBEFORE=false;
    end
    userfun = @(x) ChokeBoundModelFunctionSNOPT(x);
    
    nDim = 2;
    saveX = zeros([nDim,0]);

    dim = ones(nDim,1);
    
    if nargin == 0;
        startPt = 0*dim;
    elseif numel(startPt)==1;
        startPt = (startPt) * dim;
    end
    if nargin<4
        plotFormat = {'linewidth',2,'marker','o'};
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
    
    if isConstrAct
        plotFormat(end+1:end+2) = {'linestyle','--'};
    end
    
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
    plotPoints = @(points,form)plot(points([1:end],1),points([1:end],2),form{:});
    if isNewPlot
        h = figure;
        hold on;
        
        [x,y] = meshgrid(unique(linspace(-2,2,100)));
        z = reshape(ChokeBoundModelFunctionSNG([x(:),y(:)]),size(x));
        contour(x,y,z,30)
    end
    l = plotPoints(saveX', plotFormat);
    constrStr = 'noconstr';
    if isConstrAct
        constrStr = 'constr';
    end
    l.DisplayName = [num2str(reshape(startPt,[1,numel(startPt)])),'; ',...
        constrStr, '; f*=', num2str(fval)];
    plot(xEnd(1),xEnd(2),'Color', l.Color,'markersize',10,'marker','x', 'linewidth',2);
    ls = findobj(l.Parent.Parent,'type','line');
    legend(ls(~cellfun(@isempty,{ls.DisplayName})));
    axis([-3 3 -3 3])
    function [out, G]=ChokeBoundModelFunctionSNOPT(x)
        saveX(:,end+1)=x;
        f = @(x) 0.5*sum(x,1) + sum(x.^2,1) + (sum(x.^2,1)<1).*(10- 3*sum(x.^2,1));
        out = f(x);
        G = zeros(size(x));
        eps = 1e-8;
        deltaX = zeros(size(x));
        for ii = 1:numel(size(x))
            deltaX(ii) = eps;
            G(ii) = (f(x+deltaX) - f(x-deltaX))/(2*eps);
            deltaX(ii)=0;
        end
    end
    
    function [c, ceq, gc, gceq] = nlcon(x)
        c = -sum(x.^2)+isConstrAct*2-1;
        ceq = [];
        gc = -2.*x;
        gceq = [];
    end
    SnoptChokeSim_HASRUNBEFORE = true;
end


function [out]=ChokeBoundModelFunctionSNG(x)

        f = @(x) 0.5*sum(x,2) + sum(x.^2,2) + (sum(x.^2,2)<1).*(10- 3*sum(x.^2,2));
        out = f(x);


end





