function []=Run_SnoptChokeSim(is3D)

    pts = [ 2 2
        -2 2
        -2 -2
        2 -2
        ];
    h = figure('Name','Test Discontinuity vs constraint vs Failure on SNOPT');
    
    ax = subplot(1,4,1);
    [x,y] = meshgrid(unique(linspace(-2,2,60)));
    z = reshape(ChokeBoundModelFunctionNoGrad([x(:),y(:)]),size(x));
    surf(x,y,z)
    zlabel('$J(x)$')
    kStart = 1;
    kk = 1;
    axisNames = {'No constraint', 'Constraint $\sum x_i^2 \geq 1$', 'Failure'};
    for isConstr = [0:1,nan]
        isNewPlot = 1;
        ax(kk) = subplot(1,4,kStart+kk);
        for pt = pts'
            SnoptChokeSim(pt, isConstr,isNewPlot,is3D, ax(kk));
            isNewPlot = 0;
        end
        kk= kk+1;
    end
    for ii =1:numel(ax)
        ax(ii).Title.String = axisNames{ii};
    end
    try
        FigureTextTools(h)
        activeInterpreter = 'latex';
    catch
    end
    legs = findobj(h,'type','legend');
    [legs.Location]=deal('SouthOutside');
    ls = findobj(h,'type','line');
    newDisp = regexprep({ls.DisplayName},';.*constr.*;',':');
    [ls.DisplayName]=deal(newDisp{:});
    newDisp = regexprep({ls.DisplayName},''')','');
    [ls.DisplayName]=deal(newDisp{:});
    for ii = 1:numel(legs)
    legs(ii).Title.String = '$(\mathbf{x}_0)\rightarrow(\mathbf{x}^*): f^*$';
    legs(ii).Title.Interpreter = 'latex';
    end
    [legs.Location]=deal('SouthOutside');
    axs = findobj(h,'type','axes');
    for ii = 1:numel(axs);
        axs(ii).XLabel.String= 'X';
    	axs(ii).YLabel.String= 'Y';
        axs(ii).Title.Interpreter= 'latex';
    end

    for l = ls'
    if ~strcmp(l.LineStyle, 'none')
    l.LineWidth = 1;
    end
    end
end

function [out]=ChokeBoundModelFunctionNoGrad(x)

        f = @(x) 0.5*sum(x,2) + sum(x.^2,2) + (sum(x.^2,2)<1).*(10- 3*sum(x.^2,2));
        out = f(x);


end