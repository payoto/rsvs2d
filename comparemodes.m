
%% functions
plotPoints= @(points,str) plot(points([1:end],1),points([1:end],2));
naca4t=@(x,t)  5*t*(0.2969*sqrt(x)-0.1260*x-0.3516*x.^2+0.2843*x.^3-0.1036*x.^4);
naca4c=@(x,m,p) [m*x(x<=p)/p^2.*(2*p-x(x<=p)) , m*(1-x(x>=p))/(1-p)^2.*(1+x(x>=p)-2*p)];
%% Data extraction
pathStr1='C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\results\Optimisation\HPC\InverseDesign\NACA2212\Dir_2016-09-26T115519_bp32_2212cos\iteration_33';

pathStr2='C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\results\Optimisation\HPC\InverseDesign\NACA2212\Dir_2016-09-26T115519_bp32_2212cos\iteration_34';

[sensloop1,sensloop2]=ValidateSensivity(pathStr1,pathStr2);
%% figures

x=linspace(0,1,2000);
figure
hold on
plotPoints([[x';flip(x)'],[naca4c(x,2/100,2/10)'+naca4t(x',.12);...
    flip(naca4c(x,2/100,2/10)'-naca4t(x',0.12))]])

for ii=1:numel(sensloop1),
    for jj=1:numel(sensloop1(ii).loop),
        plotPoints(sensloop1(ii).loop(jj).subdivspline),
    end,
end


figure
hold on
plotPoints([[x';flip(x)'],[naca4c(x,2/100,2/10)'+naca4t(x',.12);...
    flip(naca4c(x,2/100,2/10)'-naca4t(x',0.12))]])

for ii=1:numel(sensloop2),
    for jj=1:numel(sensloop2(ii).loop),
        plotPoints(sensloop2(ii).loop(jj).subdivspline),
    end,
end

%% Difference figures
figure
hold on

for ii=1:numel(sensloop1)/2,
    dist=sqrt(sum((sensloop1(ii).loop.subdivspline-sensloop1(1).loop.subdivspline).^2,2));
    diff=sensloop1(ii).loop.subdivspline-sensloop1(1).loop.subdivspline;
    plotPoints([sensloop1(ii).loop.subdivspline(:,1),diff(:,2)]),
    
end