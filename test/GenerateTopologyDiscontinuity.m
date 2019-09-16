function [h,ls]=GenerateTopologyDiscontinuity(out3)
    h=[];
    ls=[];
    if nargin==0
        [h(1)]=Plot_SingleVOS();
        [h(2)]=Plot_DoubleVOS();
        [h]=Plot_DoublePatchVOS();
    elseif nargin==1
        [h,ls]=Plot_TopologyChange2Fig(out3);
        [h,ls]=Plot_TopologyPatch2Fig(out3);
    end
end

function [h]=Plot_SingleVOS()
    boxP = [1 0 0;
        2 0 0;
        2 1 0;
        1 1 0;
        1 0 0];

    V1 = [1 0; 2 0; 2 1; 1 1; 1 0];
    n = size(V1,1);
    nRep = 20;
    pts = zeros(0,2); % nRep*(n-1)+1
    for ii = 1:n-1
        pts = [pts; ones(nRep,1)*V1(ii,:)];
    end
    pts = [pts; V1(end,:)];
    theta= linspace(-pi,pi,size(pts,1))';
    C1 = [cos(theta)/2+1.5, sin(theta)/2+0.5];

    h=figure('name','SingleVOSProgression','position',[300 300 1000 600]);
    
    subplot(1, 2, 1)
    axis equal
    hold on, 
    view(-10, 20)
    plot3(boxP(:,1),boxP(:,2),boxP(:,3),'k--','linewidth',2,'DisplayName', 'VOS boundary')
    plot3(pts(:,1),pts(:,2), ones(size(pts(:,1))),'linewidth',1, 'DisplayName','$\left( 1\right)$')
    plot3(C1(:,1),C1(:,2), 0.7854*ones(size(C1(:,1))),'linewidth',1, 'DisplayName','$\left(0.79\right)$')
    plot3(1.5*ones(size(C1(:,1))), 0.5*ones(size(C1(:,1))), zeros(size(C1(:,1))),'linewidth',1, 'DisplayName','$\left(0\right)$')
    l=legend(flip(findobj(h,'type','line')), 'location', 'southeast');
    l.Title.String = 'Design variable $\left(a_1\right)$';
        
    ylabel('Y')
    zlabel('$a_1$')

    subplot(1, 2, 2)
    axis equal
    hold on, 
    view(-10, 20)

    surf([pts(:,1) C1(:,1) 1.5*ones(size(C1(:,1)))],...
        [pts(:,2) C1(:,2) 0.5*ones(size(C1(:,1)))], ...
        [ones(size(pts(:,1))) 0.7854*ones(size(C1(:,1))) 0*ones(size(C1(:,1)))],...
        'facecolor',[0.5,0.5,0.5], 'edgecolor',[0.7,0.7,0.7])
    plot3(boxP(:,1),boxP(:,2),boxP(:,3),'k--','linewidth',2,'DisplayName', 'VOS boundary')
    plot3(pts(:,1),pts(:,2), ones(size(pts(:,1))),'linewidth',3, 'DisplayName','$\left( 1\right)$')
    plot3(C1(:,1),C1(:,2), 0.7854*ones(size(C1(:,1))),'linewidth',3, 'DisplayName','$\left(0.79\right)$')
    plot3(1.5*ones(size(C1(:,1))), 0.5*ones(size(C1(:,1))), zeros(size(C1(:,1))),'linewidth',3, 'DisplayName','$\left(0\right)$')
    xlabel('X')


    try
        FigureTextTools(h,'fontSizeDelta',2)
    catch
    end
end

function [h]=Plot_DoubleVOS()
         
    boxP = [0 0 0;
        2 0 0;
        2 1 0;
        0 1 0;
        0 0 0;
        nan nan nan;
        1 0 0;
        1 1 0];
    
    V1 = [1 0; 2 0; 2 1; 1 1; 1 0];
    n = size(V1,1);
    nRep = 20;
    pts = zeros(0,2); % nRep*(n-1)+1
    for ii = 1:n-1
        pts = [pts; ones(nRep,1)*V1(ii,:)];
    end
    
    %arc from 1 1 to 1 0 with centre at [3 0.5]
    xCirc = @(y) 3 - sqrt((4.25)-(y-0.5).^2);
    yC = linspace(1,0,20)';
    pCirc = [xCirc(yC), yC];
        
    theta= [linspace(-pi,pi,size(pts-20,1))';-pi*ones(size(pCirc,1),1)];
    
    c2Size = 0.075;
    pts = [pts; pCirc];
    C1 = [cos(theta)/2+1.4, sin(theta)/2+0.5];
    C1_5 = [cos(theta)/2+1.4, sin(theta)/2+0.5];
    C1_5(theta>-pi/2 & theta<pi/2,1) = C1_5(theta>-pi/2 & theta<pi/2,1)+0.1;
    C2 = [cos(theta)*c2Size+(1-c2Size), sin(theta)*c2Size+0.5];
    
    C4 = [cos(theta)*c2Size+0.5, sin(theta)*c2Size+0.5];
    C3 = [cos(theta)*c2Size+0.5+(0.5-c2Size)/2, sin(theta)*c2Size+0.5];
    
    c2Size = 0.075/2;
    C5 = [cos(theta)*c2Size+0.5, sin(theta)*c2Size+0.5];

    h=figure('name','TwoVOSProgression','position',[300 300 1000 600]);
    subplot(1, 2, 1)
    axis equal
    hold on, 
    plot3(pts(:,1),pts(:,2), ones(size(pts(:,1))),'linewidth',1, 'DisplayName','$\left( 0.2,\ 1\right)$')
    plot3(C1_5(:,1),C1_5(:,2), 0.85*ones(size(C1_5(:,1))),'linewidth',1, 'DisplayName','$\left( 0.2,\ 0.85\right)$')
    plot3(C1(:,1),C1(:,2), 0.7854*ones(size(C1(:,1))),'linewidth',1, 'DisplayName','$\left( 0.2,\ 0.7854 \right)$')
    plot3(C2(:,1),C2(:,2), 0*ones(size(C1(:,1))),'linewidth',1, 'DisplayName','$\left( 0.2,\ +\rightarrow 0\right)$')
    plot3(C4(:,1),C4(:,2), -0*ones(size(C1(:,1))),'linewidth',1, 'DisplayName','$\left( 0.2, 0\right)$')
    plot3(C5(:,1),C5(:,2), -0.2*ones(size(C1(:,1))),'linewidth',1, 'DisplayName','$\left( 0.1,\ 0 \right)$')
    plot3(boxP(:,1),boxP(:,2),boxP(:,3),'k--','linewidth',2,'DisplayName', 'VOS boundaries')
    l=legend(flip(findobj(h,'type','line')), 'location', 'southeast');
    
    ylabel('Y')
    zlabel('pseudo VOS')

    l.Title.String = char('Design variables',' $\left(a_1,\ a_2\right)$');
    grid on
    view(-10, 20)
    subplot(1, 2, 2)
    axis equal
    hold on, 
    col1 = [0.5 0.5 0.7];
    col1Line = [0.7 0.7 0.7];
    
    col2 = [0.5 0.7 0.5];
    colpatch = [0.7 0.1 0.1];
    surf([pts(:,1)  C1_5(:,1)  C1(:,1) C2(:,1)],...
        [pts(:,2) C1_5(:,2)  C1(:,2) C2(:,2)], ...
        [ones(size(pts(:,1))) 0.85*ones(size(C1(:,1))) 0.7854*ones(size(C1(:,1))) 0*ones(size(C1(:,1))) ],...
        'facecolor',col1, 'edgecolor',col1Line,...
        'DisplayName','Change of $a_2$')
     % surf([C2(:,1),C3(:,1), C4(:,1)],...
     %    [C2(:,2), C3(:,2), C4(:,2)], ...
     %    [0*ones(size(C1(:,1))) -0*ones(size(C1(:,1))) -0*ones(size(C1(:,1)))],...
     %    'linestyle', ':', 'facecolor',colpatch, 'edgecolor',col1Line)
    surf([C4(:,1), C5(:,1)],...
        [C4(:,2), C5(:,2)], ...
        [-0*ones(size(C1(:,1))) -0.2*ones(size(C1(:,1)))],...
        'facecolor',col2, 'edgecolor',col1Line,'linestyle','none',...
        'DisplayName','Change of $a_1$')
    
    l=legend(flip(findobj(h,'type','surf')), 'location', 'southeast');
    
    plot3(pts(:,1),pts(:,2), ones(size(pts(:,1))),'linewidth',3)
    plot3(C1_5(:,1),C1_5(:,2), 0.85*ones(size(C1_5(:,1))),'linewidth',3)
    plot3(C1(:,1),C1(:,2), 0.7854*ones(size(C1(:,1))),'linewidth',3)
    plot3(C2(:,1),C2(:,2), 0*ones(size(C1(:,1))),'linewidth',3)
    plot3(C4(:,1),C4(:,2), -0.0*ones(size(C1(:,1))),'linewidth',3)
    plot3(C5(:,1),C5(:,2), -0.2*ones(size(C1(:,1))),'linewidth',3)

    plot3(boxP(:,1),boxP(:,2),boxP(:,3),'k--','linewidth',2)
    xlabel('X')
    
    grid on
    view(-10, 20)
    try
        FigureTextTools(h,'fontSizeDelta',2)
    catch
    end
end

function [h]=Plot_DoublePatchVOS()
        
    boxP = [0 0 0;
        2 0 0;
        2 1 0;
        0 1 0;
        0 0 0;
        nan nan nan;
        1 0 0;
        1 1 0];
    
    V1 = [1 0; 2 0; 2 1; 1 1; 1 0];
    n = size(V1,1);
    nRep = 20;
    pts = zeros(0,2); % nRep*(n-1)+1
    for ii = 1:n-1
        pts = [pts; ones(nRep,1)*V1(ii,:)];
    end
    
    %arc from 1 1 to 1 0 with centre at [3 0.5]
    xCirc = @(y) 3 - sqrt((4.25)-(y-0.5).^2);
    yC = linspace(1,0,20)';
    pCirc = [xCirc(yC), yC];
        
    theta= [linspace(-pi,pi,size(pts-20,1))';-pi*ones(size(pCirc,1),1)];
    
    c2Size = 0.075;
    pts = [pts; pCirc];
    C1 = [cos(theta)/2+1.4, sin(theta)/2+0.5];
    C1_5 = [cos(theta)/2+1.4, sin(theta)/2+0.5];
    C1_5(theta>-pi/2 & theta<pi/2,1) = C1_5(theta>-pi/2 & theta<pi/2,1)+0.1;
    C2 = [cos(theta)*c2Size+(1-c2Size), sin(theta)*c2Size+0.5];
    
    C4 = [cos(theta)*c2Size+0.5, sin(theta)*c2Size+0.5];
    C3 = [cos(theta)*c2Size+0.5+(0.5-c2Size)/2, sin(theta)*c2Size+0.5];
    
    c2Size = 0.075/2;
    C5 = [cos(theta)*c2Size+0.5, sin(theta)*c2Size+0.5];

    h=figure('name','TwoVOSProgression','position',[300 300 1000 600]);
    subplot(1, 2, 1)
    axis equal
    hold on, 
    plot3(boxP(:,1),boxP(:,2),boxP(:,3),'k--','linewidth',2,'DisplayName', 'VOS boundaries')
    plot3(pts(:,1),pts(:,2), ones(size(pts(:,1))),'linewidth',1, 'DisplayName','$\left( 0.2,\ 1,\ n.d. \right)$')
    plot3(C1_5(:,1),C1_5(:,2), 0.85*ones(size(C1_5(:,1))),'linewidth',1, 'DisplayName','$\left( 0.2,\ 0.85,\ n.d. \right)$')
    plot3(C1(:,1),C1(:,2), 0.7854*ones(size(C1(:,1))),'linewidth',1, 'DisplayName','$\left( 0.2,\ 0.7854,\ n.d. \right)$')
    plot3(C2(:,1),C2(:,2), 0*ones(size(C1(:,1))),'linewidth',1, 'DisplayName','$\left( 0.2,\ {+\rightarrow0},\ 1 \right)$')
    plot3(C3(:,1),C3(:,2), -0.25*ones(size(C1(:,1))),'linewidth',1, 'DisplayName','$\left( 0.1,\ 0 ,\ 0.5\right)$', 'color',[0.0784    0.1686    0.5490])
    plot3(C4(:,1),C4(:,2), -0.5*ones(size(C1(:,1))),'linewidth',1, 'DisplayName','$\left( 0.2,\ 0,\ 0 \right)$')
    plot3(C5(:,1),C5(:,2), -0.7*ones(size(C1(:,1))),'linewidth',1, 'DisplayName','$\left( 0.1,\ 0,\ n.d. \right)$')
    l=legend(flip(findobj(h,'type','line')), 'location', 'southeast');
    ylabel('Y')
    zlabel('pseudo VOS')

    l.Title.String = char('Design variables',' $\left(a_1,\ a_2,\ \sigma\right)$');
    grid on
    view(-10, 20)
    subplot(1, 2, 2)
    axis equal
    hold on, 
    col1 = [0.5 0.5 0.7];
    col1Line = [0.7 0.7 0.7];
    
    col2 = [0.5 0.7 0.5];
    colpatch = [0.7 0.1 0.1];
    surf([pts(:,1)  C1_5(:,1)  C1(:,1) C2(:,1)],...
        [pts(:,2) C1_5(:,2)  C1(:,2) C2(:,2)], ...
        [ones(size(pts(:,1))) 0.85*ones(size(C1(:,1))) 0.7854*ones(size(C1(:,1))) 0*ones(size(C1(:,1))) ],...
        'facecolor',col1, 'edgecolor',col1Line,...
        'DisplayName','Change of $a_2$')
    surf([C2(:,1),C3(:,1), C4(:,1)], ...
        [C2(:,2), C3(:,2), C4(:,2)], ...
        [0*ones(size(C1(:,1))) -0.25*ones(size(C1(:,1))) -0.5*ones(size(C1(:,1)))],...
        'linestyle', ':', 'facecolor',colpatch, 'edgecolor',col1Line,...
        'DisplayName','Change of $\sigma$')
    surf([C4(:,1), C5(:,1)], ...
        [C4(:,2), C5(:,2)], ...
        [-0.50*ones(size(C1(:,1))) -0.7*ones(size(C1(:,1)))],...
        'facecolor',col2, 'edgecolor',col1Line,'linestyle','none',...
        'DisplayName','Change of $a_1$')

    l=legend(flip(findobj(h,'type','surf')), 'location', 'southeast');
    
    plot3(pts(:,1),pts(:,2), ones(size(pts(:,1))),'linewidth',3)
    plot3(C1_5(:,1),C1_5(:,2), 0.85*ones(size(C1_5(:,1))),'linewidth',3)
    plot3(C1(:,1),C1(:,2), 0.7854*ones(size(C1(:,1))),'linewidth',3)
    plot3(C2(:,1),C2(:,2), 0*ones(size(C1(:,1))),'linewidth',3)
    plot3(C3(:,1),C3(:,2), -0.25*ones(size(C1(:,1))),'linewidth',3, 'color',[0.0784    0.1686    0.5490])
    plot3(C4(:,1),C4(:,2), -0.5*ones(size(C1(:,1))),'linewidth',3)
    plot3(C5(:,1),C5(:,2), -0.7*ones(size(C1(:,1))),'linewidth',3)

    plot3(boxP(:,1),boxP(:,2),boxP(:,3),'k--','linewidth',2)
    xlabel('X')
    grid on
    view(-10, 20)
    try
        FigureTextTools(h,'fontSizeDelta',2)
    catch
    end
end

function [h,ls]=Plot_TopologyChange2Fig(out3)
    
    
    h=figure('name','TopologySplit','position',[300 300 1000 600]);
    ax(1)=subplot(2,2,1);
    hold on;
    Plot_TopologyChange(out3, 0);
    ax(2)=subplot(2,2,2);
    hold on;
    Plot_TopologyChange(out3, 1);
    
    ax(3)=subplot(2,2,3);
    hold on;
    Plot_TopologyChange(out3, 0);
    ax(4)=subplot(2,2,4);
    hold on;
    Plot_TopologyChange(out3, 1);
    FigureTextTools(h,'fontSizeDelta',2)

    ls(1)=linkprop(ax(1:2),{'CameraPosition','CameraUpVector',...
        'PlotBoxAspectRatio','XLim','ZLim','YLim'});
    ls(2)=linkprop(ax(3:4),{'CameraPosition','CameraUpVector',...
        'PlotBoxAspectRatio','XLim','ZLim','YLim'});

end

function [h,ls]=Plot_TopologyPatch2Fig(out3)
    
    
    h=figure('name','TopologyPatch','position',[300 300 1000 600]);
    ax(1)=subplot(1,2,1);
    hold on;
    Plot_TopologyPatch(out3, 0);
    ax(2)=subplot(1,2,2);
    hold on;
    Plot_TopologyPatch(out3, 1);
    
%     ax(3)=subplot(2,2,3);
%     hold on;
%     Plot_TopologyChange(out3, 0);
%     ax(4)=subplot(2,2,4);
%     hold on;
%     Plot_TopologyChange(out3, 1);
     FigureTextTools(h,'fontSizeDelta',2)

    ls(1)=linkprop(ax(1:2),{'CameraPosition','CameraUpVector',...
        'PlotBoxAspectRatio','XLim','ZLim','YLim'});
%     ls(2)=linkprop(ax(3:4),{'CameraPosition','CameraUpVector',...
%         'PlotBoxAspectRatio','XLim','ZLim','YLim'});
    
    lights = findobj(h,'type','light');
    delete(lights)

end

function []=Plot_TopologyChange(out3, withSurf)
    nPts = 50;
    boxP = [0 0 0;
        3 0 0;
        3 2 0;
        0 2 0;
        0 0 0;
        nan nan nan;
        1 0 0;
        1 2 0;
        nan nan nan;
        2 0 0;
        2 2 0];
    nStart = 3;
    
    plot3(boxP(:,1),boxP(:,2),boxP(:,3),'k--','linewidth',2,'DisplayName', 'VOS boundaries')
    % Split around LE TE
    for ii = 1:numel(out3)
        for jj = 1:numel(out3(ii).plotdat)
            isCCW = CCWLoop(out3(ii).plotdat(jj).ctrl);
            if ~isCCW
                out3(ii).plotdat(jj).ctrl = flip(out3(ii).plotdat(jj).ctrl);
            end
            [~,iMin]=min(out3(ii).plotdat(jj).ctrl(:,1)+abs(out3(ii).plotdat(jj).ctrl(:,2)-1));
            if ii < 10
                out3(ii).plotdat(jj).ctrl = out3(ii).plotdat(jj).ctrl([iMin:end-1,1:iMin-1],:);
            else
                out3(ii).plotdat(jj).ctrl = out3(ii).plotdat(jj).ctrl([iMin:end,1:iMin-1],:);
            end
            
            out3(ii).plotdat(jj).pts=RemoveIdenticalConsecutivePoints(...
                out3(ii).plotdat(jj).pts);
            if ~CCWLoop(out3(ii).plotdat(jj).pts)
                out3(ii).plotdat(jj).pts = flip(out3(ii).plotdat(jj).pts);
            end
            [~,iMin]=min(out3(ii).plotdat(jj).pts(:,1)+abs(out3(ii).plotdat(jj).pts(:,2)-1));
            out3(ii).plotdat(jj).pts = out3(ii).plotdat(jj).pts([iMin:end,1:iMin-1],:);
            
        end
    end
    ii = 10; jj = 1;
    out3(ii).plotdat(jj).pts = flip(out3(ii).plotdat(jj).pts);
    % Cage
    DatCage = zeros(3,1);
    jj=1;
    ll = 1;
    for kk = 1:2:size(out3(1).plotdat(jj).ctrl,1)
        for ii = nStart:9
            DatCage(ll,1:2)=out3(ii).plotdat(jj).ctrl(kk,:);
            DatCage(ll,3)=out3(ii).unstructured.cell.fill(8);
            ll = ll+1;
        end
        DatCage(ll,1:3) = [nan nan nan];
        ll = ll+1;
    end
    
    for ii = 1:11
        for jj = 1:numel(out3(ii).plotdat)
            out3(ii).plotdat(jj).pts=SplitAzimuthPatch(out3(ii).plotdat(jj).pts,nPts);
        end
    end
    
    jj=1;
    
    surfCageX = zeros(size(out3(ii).plotdat(jj).pts([1:end,1],1),1),0);
    surfCageY = surfCageX;
    surfCageZ = surfCageX;
    maxN = numel(out3(1).plotdat(1).pts([1:end],2));
    
    for ii = nStart:9
        surfCageX=[surfCageX,out3(ii).plotdat(jj).pts([1:maxN,1],1)];
        surfCageY=[surfCageY,out3(ii).plotdat(jj).pts([1:maxN,1],2)];
        surfCageZ=[surfCageZ,ones(size(out3(ii).plotdat(jj).pts([1:maxN,1],2)))*out3(ii).unstructured.cell.fill(8)];

    end
    if withSurf
    s(1) = surf(surfCageX,surfCageY,surfCageZ,'linestyle','none',...
        'facecolor',[0.75 0.75 0.75],...
        'DisplayName','Smooth $a_2$ variation');
    end
    plot3(DatCage(:,1),DatCage(:,2),DatCage(:,3),'k-','linewidth',1,'DisplayName', '')
    for ii = nStart:9
        jj = 1;
        l(ii) = plot3(out3(ii).plotdat(jj).pts(:,1),...
            out3(ii).plotdat(jj).pts(:,2),...
            out3(ii).unstructured.cell.fill(8)*ones(size(out3(ii).plotdat(jj).pts(:,2))),'-',...
            'DisplayName', sprintf('($0.1,\\ %.2f,\\ 0.1$)',out3(ii).unstructured.cell.fill(8)),...
            'linewidth',2);
%         l2(ii) = plot3(out3(ii).plotdat(jj).ctrl(1:end,1),...
%             out3(ii).plotdat(jj).ctrl(1:end,2),...
%             out3(ii).unstructured.cell.fill(8)*ones(size(out3(ii).plotdat(jj).ctrl(1:end,2))),'s',...
%             'DisplayName', int2str(ii));
%         l2(ii).Color = l(ii).Color;
        
    end
    
    % Make the 12th
    [out3]=ReplicateAndScale(out3,11,12,0.7,[-0.12 0.12],0.18,-0.1);
    [out3]=ReplicateAndScale(out3,11,11,1.1,[-0.06 0.06],0.28,0);
    out3(11).unstructured.cell.fill(7) = out3(12).unstructured.cell.fill(7);
    out3(11).unstructured.cell.fill(8) = out3(10).unstructured.cell.fill(8)*0.95;

    n = 0;
    for jj = 1:numel(out3(10).plotdat)
        n = n + size(out3(10).plotdat(jj).pts([1:end,1],1),1)+1;
    end
    
    surfCageX = zeros(n,0);
    surfCageY = surfCageX;
    surfCageZ = surfCageX;
    
    for ii = 10:11
        currPts = zeros(0,2);
        for jj = 1:numel(out3(ii).plotdat)
            currPts = [currPts; out3(ii).plotdat(jj).pts([1:end,1],:); [nan nan]];
        end
        surfCageX=[surfCageX, currPts(:,1)];
        surfCageY=[surfCageY, currPts(:,2)];
        surfCageZ=[surfCageZ, ones(size(currPts(:,1)))*out3(ii).unstructured.cell.fill(8)];
    end
    if withSurf
    s(2) = surf(surfCageX,surfCageY,surfCageZ,'linestyle','none',...
        'facecolor',[0.8 0.4 0.4],...
        'DisplayName','Discontinuity in $a_2$');
    end
    surfCageX = zeros(n,0);
    surfCageY = surfCageX;
    surfCageZ = surfCageX;
    
    
    % out3(12) = out3(11);
    % out3(12).unstructured.cell.fill(8) = -0.2;
    % d = [-0.12 0.12];
    % for jj = 1:numel(out3(11).plotdat)
    %     repThisMat = @(vec) repmat(vec,[size(out3(11).plotdat(jj).pts,1),1]);
    %     deltaCentre = repThisMat(mean(out3(11).plotdat(jj).pts));
    %     out3(12).plotdat(jj).pts = (out3(12).plotdat(jj).pts-deltaCentre)*0.7...
    %         +repThisMat(mean(out3(11).plotdat(jj).pts)+[d(jj) 0]);
        
    %     repThisMat = @(vec) repmat(vec,[size(out3(11).plotdat(jj).ctrl,1),1]);
    %     deltaCentre = repThisMat(mean(out3(11).plotdat(jj).pts));
    %     out3(12).plotdat(jj).ctrl = (out3(12).plotdat(jj).ctrl-deltaCentre)*0.7...
    %         +repThisMat(mean(out3(11).plotdat(jj).pts)+[d(jj) 0 ]);
    % end
    % delExtra = 0.18;
    % out3(12).plotdat(2).ctrl(2,:) = [1 1-delExtra ];
    % out3(12).plotdat(2).ctrl(4,:) = [1 1+delExtra];
    % out3(12).plotdat(1).ctrl(2,:) = [2 1-delExtra ];
    % out3(12).plotdat(1).ctrl(4,:) = [2 1+delExtra];
    % out3(12).unstructured.cell.fill(7) = out3(12).unstructured.cell.fill(7)*0.7^2;
    
    % make the datcage
    
    for ii = 11:12
        currPts = zeros(0,2);
        for jj = 1:numel(out3(11).plotdat)
            currPts = [currPts; out3(ii).plotdat(jj).pts([1:end,1],:); [nan nan]];
        end
        surfCageX=[surfCageX, currPts(:,1)];
        surfCageY=[surfCageY, currPts(:,2)];
        surfCageZ=[surfCageZ, ones(size(currPts(:,1)))*out3(ii).unstructured.cell.fill(8)];
    end
    if withSurf
        s(3) = surf(surfCageX,surfCageY,surfCageZ,'linestyle','none',...
            'facecolor',[0.4 0.4 0.8],...
            'DisplayName','Smooth $a_2$ variation after cut');
    end
    surfCageX = zeros(n,0);
    surfCageY = surfCageX;
    surfCageZ = surfCageX;
    for ii = 12:numel(out3)
        currPts = zeros(0,2);
        for jj = 1:numel(out3(11).plotdat)
            currPts = [currPts; out3(ii).plotdat(jj).pts([1:end,1],:); [nan nan]];
        end
        surfCageX=[surfCageX, currPts(:,1)];
        surfCageY=[surfCageY, currPts(:,2)];
        surfCageZ=[surfCageZ, ones(size(currPts(:,1)))*out3(ii).unstructured.cell.fill(8)];
    end
    if withSurf
    s(4) = surf(surfCageX,surfCageY,surfCageZ,'linestyle','none',...
        'facecolor',[0.4 0.8 0.4],...
        'DisplayName','Variation in $a_1$ and $a_3$');
    end
    DatCage = zeros(3,1);
    ll = 1;
    for jj = 1:numel(out3(10).plotdat)
        for kk = 1:size(out3(10).plotdat(jj).ctrl,1)
            for ii = 10:numel(out3)
                DatCage(ll,1:2)=out3(ii).plotdat(jj).ctrl(kk,:);
                DatCage(ll,3)=out3(ii).unstructured.cell.fill(8);
                ll = ll+1;
            end
            DatCage(ll,1:3) = [nan nan nan];
            ll = ll+1;
        end
    end
    plot3(DatCage(:,1),DatCage(:,2),DatCage(:,3),'k-','linewidth',1,'DisplayName', '')

    
    for ii = 11:numel(out3)
        for jj = 1:size(out3(ii).loopNurbs,2)
            if jj==1
                l(ii) = plot3(out3(ii).plotdat(jj).pts(:,1),...
                    out3(ii).plotdat(jj).pts(:,2),...
                    out3(ii).unstructured.cell.fill(8)*ones(size(out3(ii).plotdat(jj).pts(:,2))),...
                    'DisplayName', sprintf('($%.2f,\\ %.2f,\\ %.2f$)',out3(ii).unstructured.cell.fill(7),...
                    out3(ii).unstructured.cell.fill(8),...
                    out3(ii).unstructured.cell.fill(7)),...
                    'linewidth',2);
            else
                plot3(out3(ii).plotdat(jj).pts(:,1),...
                    out3(ii).plotdat(jj).pts(:,2),...
                    out3(ii).unstructured.cell.fill(8)*ones(size(out3(ii).plotdat(jj).pts(:,2))),...
                    'DisplayName', sprintf('(0.1,\\ $%.2f,\\ 0.1$)',out3(ii).unstructured.cell.fill(8)),...
                    'linewidth',2, 'Color', l(ii).Color);
            end
%             l2(ii) = plot3(out3(ii).plotdat(jj).ctrl(1:2:end,1),...
%                 out3(ii).plotdat(jj).ctrl(1:2:end,2),...
%                 out3(ii).unstructured.cell.fill(8)*ones(size(out3(ii).plotdat(jj).ctrl(1:2:end,2))),'s',...
%                 'DisplayName', int2str(ii));
%             l2(ii).Color = l(ii).Color;
        end
    end

    if withSurf

        legO=legend([s]);
%         legO.Title.String=char('V', '($a_1,\ a_2,\ a_3$)');
    else
        legO=legend([l([nStart:9,11:end])]);
        legO.Title.String=char('VOS values', '($a_1,\ a_2,\ a_3$)');
    end
    lObj = light();
%     axis equal
end 

function []=Plot_TopologyPatch(out3, withSurf)
    nPts = 50;
    boxP = [0 0 0;
        3 0 0;
        3 2 0;
        0 2 0;
        0 0 0;
        nan nan nan;
        1 0 0;
        1 2 0;
        nan nan nan;
        2 0 0;
        2 2 0];
    nStart = 3;
    
    plot3(boxP(:,1),boxP(:,2),boxP(:,3),'k--','linewidth',2,'DisplayName', 'VOS boundaries')
    % Split around LE TE
    for ii = 1:numel(out3)
        for jj = 1:numel(out3(ii).plotdat)
            isCCW = CCWLoop(out3(ii).plotdat(jj).ctrl);
            if ~isCCW
                out3(ii).plotdat(jj).ctrl = flip(out3(ii).plotdat(jj).ctrl);
            end
            [~,iMin]=min(out3(ii).plotdat(jj).ctrl(:,1)+abs(out3(ii).plotdat(jj).ctrl(:,2)-1));
            if ii < 10
                out3(ii).plotdat(jj).ctrl = out3(ii).plotdat(jj).ctrl([iMin:end-1,1:iMin-1],:);
            else
                out3(ii).plotdat(jj).ctrl = out3(ii).plotdat(jj).ctrl([iMin:end,1:iMin-1],:);
            end
            
            out3(ii).plotdat(jj).pts=RemoveIdenticalConsecutivePoints(...
                out3(ii).plotdat(jj).pts);
            if ~CCWLoop(out3(ii).plotdat(jj).pts)
                out3(ii).plotdat(jj).pts = flip(out3(ii).plotdat(jj).pts);
            end
            [~,iMin]=min(out3(ii).plotdat(jj).pts(:,1)+abs(out3(ii).plotdat(jj).pts(:,2)-1));
            out3(ii).plotdat(jj).pts = out3(ii).plotdat(jj).pts([iMin:end,1:iMin-1],:);
            
        end
    end
    ii = 10; jj = 1;
    out3(ii).plotdat(jj).pts = flip(out3(ii).plotdat(jj).pts);
    % Cage
    DatCage = zeros(3,1);
    jj=1;
    ll = 1;
    for kk = 1:2:size(out3(1).plotdat(jj).ctrl,1)
        for ii = nStart:9
            DatCage(ll,1:2)=out3(ii).plotdat(jj).ctrl(kk,:);
            DatCage(ll,3)=out3(ii).unstructured.cell.fill(8);
            ll = ll+1;
        end
        DatCage(ll,1:3) = [nan nan nan];
        ll = ll+1;
    end
    
    for ii = 1:11
        for jj = 1:numel(out3(ii).plotdat)
            out3(ii).plotdat(jj).pts=SplitAzimuthPatch(out3(ii).plotdat(jj).pts,nPts);
        end
    end
    
    jj=1;
    
    surfCageX = zeros(size(out3(ii).plotdat(jj).pts([1:end,1],1),1),0);
    surfCageY = surfCageX;
    surfCageZ = surfCageX;
    maxN = numel(out3(1).plotdat(1).pts([1:end],2));
    
    for ii = nStart:9
        surfCageX=[surfCageX,out3(ii).plotdat(jj).pts([1:maxN,1],1)];
        surfCageY=[surfCageY,out3(ii).plotdat(jj).pts([1:maxN,1],2)];
        surfCageZ=[surfCageZ,ones(size(out3(ii).plotdat(jj).pts([1:maxN,1],2)))*out3(ii).unstructured.cell.fill(8)];

    end
    if withSurf
    s(1) = surf(surfCageX,surfCageY,surfCageZ,'linestyle','none',...
        'facecolor',[0.75 0.75 0.75],...
        'DisplayName','Smooth $a_2$ variation');
    end
    plot3(DatCage(:,1),DatCage(:,2),DatCage(:,3),'k-','linewidth',1,'DisplayName', '')
    for ii = nStart:9
        jj = 1;
        l(ii) = plot3(out3(ii).plotdat(jj).pts(:,1),...
            out3(ii).plotdat(jj).pts(:,2),...
            out3(ii).unstructured.cell.fill(8)*ones(size(out3(ii).plotdat(jj).pts(:,2))),'-',...
            'DisplayName', sprintf('($0.1,\\ %.2f,\\ 0.1$,\\ n.d.)',out3(ii).unstructured.cell.fill(8)),...
            'linewidth',2);
%         l2(ii) = plot3(out3(ii).plotdat(jj).ctrl(1:end,1),...
%             out3(ii).plotdat(jj).ctrl(1:end,2),...
%             out3(ii).unstructured.cell.fill(8)*ones(size(out3(ii).plotdat(jj).ctrl(1:end,2))),'s',...
%             'DisplayName', int2str(ii));
%         l2(ii).Color = l(ii).Color;
        
    end
    
    % Make the 12th
    [out3]=ReplicateAndScale(out3,11,12,0.7,[-0.12 0.12],0.18,-0.5);
    [out3]=ReplicateAndScale(out3,11,11,1.1,[-0.06 0.06],0.28,-0.3);
    out3(11).unstructured.cell.fill(7) = out3(12).unstructured.cell.fill(7);
    out3(11).unstructured.cell.fill(8) = out3(10).unstructured.cell.fill(8)-0.3;
    out3(12).unstructured.cell.fill(8) = out3(12).unstructured.cell.fill(8)-0.3;

    n = 0;
    for jj = 1:numel(out3(10).plotdat)
        n = n + size(out3(10).plotdat(jj).pts([1:end,1],1),1)+1;
    end
    
    surfCageX = zeros(n,0);
    surfCageY = surfCageX;
    surfCageZ = surfCageX;
    
    for ii = 10:11
        currPts = zeros(0,2);
        for jj = 1:numel(out3(ii).plotdat)
            currPts = [currPts; out3(ii).plotdat(jj).pts([1:end,1],:); [nan nan]];
        end
        surfCageX=[surfCageX, currPts(:,1)];
        surfCageY=[surfCageY, currPts(:,2)];
        surfCageZ=[surfCageZ, ones(size(currPts(:,1)))*out3(ii).unstructured.cell.fill(8)];
    end
    if withSurf
    s(2) = surf(surfCageX,surfCageY,surfCageZ,'linestyle','none',...
        'facecolor',[0.8 0.4 0.4],...
        'DisplayName','Smooth transition via $\sigma$');
    end
    surfCageX = zeros(n,0);
    surfCageY = surfCageX;
    surfCageZ = surfCageX;
    
    
    % out3(12) = out3(11);
    % out3(12).unstructured.cell.fill(8) = -0.2;
    % d = [-0.12 0.12];
    % for jj = 1:numel(out3(11).plotdat)
    %     repThisMat = @(vec) repmat(vec,[size(out3(11).plotdat(jj).pts,1),1]);
    %     deltaCentre = repThisMat(mean(out3(11).plotdat(jj).pts));
    %     out3(12).plotdat(jj).pts = (out3(12).plotdat(jj).pts-deltaCentre)*0.7...
    %         +repThisMat(mean(out3(11).plotdat(jj).pts)+[d(jj) 0]);
        
    %     repThisMat = @(vec) repmat(vec,[size(out3(11).plotdat(jj).ctrl,1),1]);
    %     deltaCentre = repThisMat(mean(out3(11).plotdat(jj).pts));
    %     out3(12).plotdat(jj).ctrl = (out3(12).plotdat(jj).ctrl-deltaCentre)*0.7...
    %         +repThisMat(mean(out3(11).plotdat(jj).pts)+[d(jj) 0 ]);
    % end
    % delExtra = 0.18;
    % out3(12).plotdat(2).ctrl(2,:) = [1 1-delExtra ];
    % out3(12).plotdat(2).ctrl(4,:) = [1 1+delExtra];
    % out3(12).plotdat(1).ctrl(2,:) = [2 1-delExtra ];
    % out3(12).plotdat(1).ctrl(4,:) = [2 1+delExtra];
    % out3(12).unstructured.cell.fill(7) = out3(12).unstructured.cell.fill(7)*0.7^2;
    
    % make the datcage
    
    for ii = 11:12
        currPts = zeros(0,2);
        for jj = 1:numel(out3(11).plotdat)
            currPts = [currPts; out3(ii).plotdat(jj).pts([1:end,1],:); [nan nan]];
        end
        surfCageX=[surfCageX, currPts(:,1)];
        surfCageY=[surfCageY, currPts(:,2)];
        surfCageZ=[surfCageZ, ones(size(currPts(:,1)))*out3(ii).unstructured.cell.fill(8)];
    end
    if withSurf
        s(3) = surf(surfCageX,surfCageY,surfCageZ,'linestyle','none',...
            'facecolor',[0.4 0.4 0.8],...
            'DisplayName','Smooth $a_2$ variation after cut');
    end
    surfCageX = zeros(n,0);
    surfCageY = surfCageX;
    surfCageZ = surfCageX;
    for ii = 12:numel(out3)
        currPts = zeros(0,2);
        for jj = 1:numel(out3(11).plotdat)
            currPts = [currPts; out3(ii).plotdat(jj).pts([1:end,1],:); [nan nan]];
        end
        surfCageX=[surfCageX, currPts(:,1)];
        surfCageY=[surfCageY, currPts(:,2)];
        surfCageZ=[surfCageZ, ones(size(currPts(:,1)))*out3(ii).unstructured.cell.fill(8)];
    end
    if withSurf
    s(4) = surf(surfCageX,surfCageY,surfCageZ,'linestyle','none',...
        'facecolor',[0.4 0.8 0.4],...
        'DisplayName','Variation in $a_1$ and $a_3$');
    end
    DatCage = zeros(3,1);
    ll = 1;
    for jj = 1:numel(out3(10).plotdat)
        for kk = 1:size(out3(10).plotdat(jj).ctrl,1)
            for ii = 10:numel(out3)
                DatCage(ll,1:2)=out3(ii).plotdat(jj).ctrl(kk,:);
                DatCage(ll,3)=out3(ii).unstructured.cell.fill(8);
                ll = ll+1;
            end
            DatCage(ll,1:3) = [nan nan nan];
            ll = ll+1;
        end
    end
    plot3(DatCage(:,1),DatCage(:,2),DatCage(:,3),'k-','linewidth',1,'DisplayName', '')

    
    for ii = 11:numel(out3)
        for jj = 1:size(out3(ii).loopNurbs,2)
            if jj==1
                l(ii) = plot3(out3(ii).plotdat(jj).pts(:,1),...
                    out3(ii).plotdat(jj).pts(:,2),...
                    out3(ii).unstructured.cell.fill(8)*ones(size(out3(ii).plotdat(jj).pts(:,2))),...
                    'DisplayName', sprintf('($%.2f,\\ %.2f,\\ %.2f,\\ n.d.$)',out3(ii).unstructured.cell.fill(7),...
                    out3(ii).unstructured.cell.fill(8),...
                    out3(ii).unstructured.cell.fill(7)),...
                    'linewidth',2);
            else
                plot3(out3(ii).plotdat(jj).pts(:,1),...
                    out3(ii).plotdat(jj).pts(:,2),...
                    out3(ii).unstructured.cell.fill(8)*ones(size(out3(ii).plotdat(jj).pts(:,2))),...
                    'DisplayName', sprintf('(0.1,\\ $%.2f,\\ 0.1,\\ 0$)',out3(ii).unstructured.cell.fill(8)),...
                    'linewidth',2, 'Color', l(ii).Color);
            end
%             l2(ii) = plot3(out3(ii).plotdat(jj).ctrl(1:2:end,1),...
%                 out3(ii).plotdat(jj).ctrl(1:2:end,2),...
%                 out3(ii).unstructured.cell.fill(8)*ones(size(out3(ii).plotdat(jj).ctrl(1:2:end,2))),'s',...
%                 'DisplayName', int2str(ii));
%             l2(ii).Color = l(ii).Color;
        end
    end

    if withSurf

        legO=legend([s]);
%         legO.Title.String=char('V', '($a_1,\ a_2,\ a_3$)');
    else
        legO=legend([l([nStart:9,11:end])]);
        legO.Title.String=char('VOS values', '($a_1,\ a_2,\ a_3,\ \sigma$)');
    end
    lObj = light();
%     axis equal
end 

function [pts]=SplitAzimuthPatch(pts,nPts)
    
    [vMin, iMin]=min(pts(:,1));
    [vMax, iMax]=max(pts(:,1));
    
    if iMin>iMax
        iMin=1;
    end
    
    p=(vMax+vMin)/2;
    indUp = iMin:iMax;
    [~,iUp] = min(abs(pts(indUp,1)-p));
    indLow = [iMax:numel(pts(:,1)), 1:iMin];
    [~,iLow] = min(abs(pts(indLow,1)-p));
    
    iUp = indUp(iUp);
    iLow = indLow(iLow);
    
    c = mod(round([linspace(iMin, iUp, nPts),linspace(iUp, iMax, nPts),...
        linspace(iMax, iLow, nPts),linspace(iLow, numel(pts(:,1)+iMin), nPts)])-1,...
        numel(pts(:,1)))+1;
    pts = pts(c,:);
end

function [out3]=ReplicateAndScale(out3,iRoot,iNew,scaleNew,shiftNew, ctrlShift, vosOffset)
    
    out3 = [out3(1:iNew-1),out3(iRoot),out3(iNew:end)];
    
    out3(iNew).unstructured.cell.fill(8) = out3(iNew).unstructured.cell.fill(8) + vosOffset;
    d = [-0.12 0.12];
    d = shiftNew;
    for jj = 1:numel(out3(iRoot).plotdat)
        repThisMat = @(vec) repmat(vec,[size(out3(iRoot).plotdat(jj).pts,1),1]);
        deltaCentre = repThisMat(mean(out3(iRoot).plotdat(jj).pts));
        out3(iNew).plotdat(jj).pts = (out3(iNew).plotdat(jj).pts-deltaCentre)*scaleNew...
            +repThisMat(mean(out3(iRoot).plotdat(jj).pts)+[d(jj) 0]);
        
        repThisMat = @(vec) repmat(vec,[size(out3(iRoot).plotdat(jj).ctrl,1),1]);
        deltaCentre = repThisMat(mean(out3(iRoot).plotdat(jj).pts));
        out3(iNew).plotdat(jj).ctrl = (out3(iNew).plotdat(jj).ctrl-deltaCentre)*scaleNew...
            +repThisMat(mean(out3(iRoot).plotdat(jj).pts)+[d(jj) 0 ]);
    end
    delExtra = ctrlShift;
    out3(iNew).plotdat(2).ctrl(2,:) = [1 1-delExtra ];
    out3(iNew).plotdat(2).ctrl(4,:) = [1 1+delExtra];
    out3(iNew).plotdat(1).ctrl(2,:) = [2 1-delExtra ];
    out3(iNew).plotdat(1).ctrl(4,:) = [2 1+delExtra];
    out3(iNew).unstructured.cell.fill(7) = out3(12).unstructured.cell.fill(7)*scaleNew^2;
end


