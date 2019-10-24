
function [saveX, discontStructAll] = SnoptChokeSim_solution(startPt, isConstrAct, isNewPlot, is3D, ax, plotFormat)
    global SnoptChokeSim_HASRUNBEFORE;
    addpath ../../ASO_LK/SRC/matlab-snopt/
    addpath ../../ASO_LK/SRC/
    if isempty(SnoptChokeSim_HASRUNBEFORE)
        SnoptChokeSim_HASRUNBEFORE=false;
    end
    % discontinuityMulti = @(x) ((tanh(100*(1-sum(x.^2,1)))+1)/2);
	discontinuityMulti = @(x) (sum(x.^2,1)<1);
    userfun = @(x) ChokeBoundModelFunctionSNOPT(x);
    
    nDim = 2;
    discontStruct.nConstr = 100;
    discontStruct.sigma= 0.5;
    discontStruct.xd = zeros(nDim,0);
    discontStruct.Dxd = discontStruct.xd ;
    discontStructAll = discontStruct;
    discontStructAll.nConstr=0;
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
    end
    ax.NextPlot = 'add';
    if ~is3D
%         l = plotPoints(ax,saveLSX', {'k--'});
        l = plotPoints(ax,saveX', plotFormat);
        plotPoints(ax,saveLSX', {plotFormat{:},'marker','s','linestyle',...
            'none','Color',l.Color});
    else
%         l = plotPoints3(ax,saveLSX', {'k--'});
        l = plotPoints3(ax,saveX', plotFormat);
    end
    if ~is3D
        plot(ax,xEnd(1),xEnd(2),'Color', l.Color,'markersize',10,'marker','x', 'linewidth',2);
    else
        plot3(ax,xEnd(1),xEnd(2),ax.ZLim(1),'Color', l.Color,'markersize',10,'marker','x', 'linewidth',2);

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
    for mm=1:discontStructAll.nConstr
        PlotConstraint(ax, discontStructAll.xd(:,mm), discontStructAll.Dxd(:,mm), ...
            discontStructAll.sigma, 'k-');
    end
    if isConstrAct
        PlotConstraint(ax, [1 0]', [-1 0]', 1, 'r-');
    end
    if ~is3D
        axis([-2 2 -2 2])
    else
        box = axis;
        box(1:4) = [-2 2 -2 2];
        axis(box)
    end

    function [out, G]=ChokeBoundModelFunctionSNOPT(x)
        
        f = @(x) 0.5*sum(x,1) + sum(x.^4+x.^3,1) + discontinuityMulti(x).*(10- 3*sum(x.^2,1));
        
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
        
        if numel(saveX)==0 || out<=saveX(end,end)
            saveX(:,end+1)=[x;out];
        else
            saveLSX(:,end+1) = [x;out];
        end
        if size(saveX,2)>=2
            if discontStruct.nConstr>0 && DetectDiscontinuity(saveX(1:end-1,end),saveX(end,end),...
                    saveX(1:end-1,end-1),saveX(end,end-1),...
                    x,out) 
                newDirection = (x-saveX(1:end-1,end))/sqrt(sum((x-saveX(1:end-1,end)).^2));
                if numel(discontStruct.xd)>0 && sum((saveX(1:end-1,end)-discontStruct.xd(:,end)).^2)<1e-5
                    discontStruct.Dxd(:,end) = (newDirection+discontStruct.Dxd(:,end))/2;
                else
                    discontStruct.xd(:,end+1) = saveX(1:end-1,end);
                    discontStruct.Dxd(:,end+1) = newDirection;
                end
                if size(discontStruct.xd,2)>discontStruct.nConstr
                    discontStruct.xd=discontStruct.xd(:,end-discontStruct.nConstr+1:end);
                    discontStruct.Dxd=discontStruct.Dxd(:,end-discontStruct.nConstr+1:end);
                end
                discontStructAll.xd(:,end+1)=discontStruct.xd(:,end);
                discontStructAll.Dxd(:,end+1)=discontStruct.Dxd(:,end);
                discontStructAll.nConstr=size(discontStructAll.xd,2);
                dbg2= true;
                if dbg2
                    PlotConstraint(gca, discontStructAll.xd(:,end), discontStructAll.Dxd(:,end), ...
                        discontStructAll.sigma, 'k-');
                    plot(saveX(1,:),saveX(2,:),'k-s','linewidth',2)
                    disp('in dbg')
                end
            end
        end
        
    end
    
    function [cAll, ceq, gcAll, gceq] = nlcon(x)
        ceq=[];
        gceq=[];
        [c, ceq, gc, gceq] = UnitCircularConstraint(x);
        %[c, gc] = GeneralCircularConstraint(x, [1 0]', [-1 0]', 1);
        cAll = zeros(1,discontStruct.nConstr+1);
        gcAll = zeros(numel(x),discontStruct.nConstr+1);
        cAll(1,1) = c;
        gcAll(:,1) = gc;
        for ii = 1:size(discontStruct.xd,2)
            [cAll(1,ii+1), gcAll(:,ii+1)]=GeneralCircularConstraint(x,...
                discontStruct.xd(:,ii), discontStruct.Dxd(:,ii), discontStruct.sigma);
        end
        for ii = 1+size(discontStruct.xd,2):discontStruct.nConstr;
            cAll(1,ii+1)=-10;
            gcAll(:,ii+1)=1;
        end
        cAll = cAll';
%         gcAll = gcAll';
    end
    function [c, ceq, gc, gceq] = UnitCircularConstraint(x)
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

            f = @(x) 0.5*sum(x,2) + sum(x.^4+x.^3,2) + discontinuityMulti(x')'.*(10- 3*sum(x.^2,2));
            
            out = f(x);
            if isnan(isConstrAct)
                out((sum(x.^2,2)<1))=nan;
            end

    end
    SnoptChokeSim_HASRUNBEFORE = true;
    dbg=false;
    
    if dbg
        [x,y] = meshgrid(unique(linspace(-2,2,100)));
        for mm=1:discontStructAll.nConstr
            Z=GeneralCircularConstraintMatrix(x,y, ...
                discontStructAll.xd(:,mm), discontStructAll.Dxd(:,mm), ...
                discontStructAll.sigma);
            Z(Z<=0)=nan;
            surf(x,y,Z)
        end 
    end
end


function Z=GeneralCircularConstraintMatrix(X,Y, xd, Dxd, sigma)
    
    Z=zeros(size(X));
    for ii =1:size(X,1)
        for jj=1:size(X,2)
            [Z(ii,jj)] = GeneralCircularConstraint([X(ii,jj);Y(ii,jj)],...
                xd, Dxd, sigma);
        end
    end
end

function [c, gc] = GeneralCircularConstraint(x, xd, Dxd, sigma)

    c = sigma.^2-sum((x-(xd+sigma*Dxd)).^2);

    gc = -2.*(x-(xd+sigma*Dxd));
end

function [isDiscontinuous] = DetectDiscontinuity(currX,currF,prevX, prevF,...
        nextX,nextF)
    % discontinuity is detected if a gradient on 1 side is 100 time larger
    % on a gap of a similar size
    
    prevDist =sqrt(sum((currX-prevX).^2));
    newDist =sqrt(sum((currX-nextX).^2));
    if newDist==0 || prevDist==0
        isDiscontinuous = false;
        return
    end
    deltaRef = abs(currF-prevF)/prevDist;
    newRef = abs(currF-nextF)/newDist;
    multiplier = 10 * max((newDist/prevDist),1);
    
    isDiscontinuous = newRef>(deltaRef*multiplier);
end

function PlotConstraint(ax, xd, Dxd, sigma, varargin)
    
    x=cos(linspace(0,2*pi));
    y=sin(linspace(0,2*pi));
    
    centre=xd+sigma*Dxd;
    plot(ax,x*sigma+centre(1),y*sigma+centre(2), varargin{:})
end





