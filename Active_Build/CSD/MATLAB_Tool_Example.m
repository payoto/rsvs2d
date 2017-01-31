
%% Create the PDE Model
% Create a pde entity for a PDE with two dependent variables
numberOfPDE = 2;
pdem = createpde(numberOfPDE);

%% Problem Parameters
E = 1.0e6; % modulus of elasticity
nu = .3; % Poisson's ratio
thick = .1; % plate thickness
len = 1.0; % side length for the square plate
hmax = len/20; % mesh size parameter
D = E*thick^3/(12*(1 - nu^2));
pres = 2; % external pressure

%% Geometry and Mesh
gdm = [3 4 0 len len 0 0 0 len len]';
g = decsg(gdm, 'S1', ('S1')');
g=[ 2 ,2 ,2 ,2 ;  % Number of elements
    5 ,10,10,0 ; % first coordinate of first vertex
    10,10,0 ,5 ; % first coordinate of second vertex
    5 ,0 ,10,10; % second coordinate of first vertex
    0 ,10,10,5 ; % second coordinate of second vertex
    1 ,1 ,1 ,1 ; % zone to the left of element (1 = inside)
    0 ,0 ,0 ,0 ]; % zone to the right of element (0=outside)
%% Create a geometry entity
geometryFromEdges(pdem,g);


% Plot the geometry and display the edge labels for use in the boundary
% condition definition.
figure;
pdegplot(pdem, 'edgeLabels', 'on');
axis equal
title 'Geometry With Edge Labels Displayed';

generateMesh(pdem, 'Hmax', hmax);
%% Boundary Conditions
k = 1e7; % spring stiffness
% Define distributed springs on all four edges
bOuter = applyBoundaryCondition(pdem,'Edge',(1:4), 'g', [0 0], 'q', [0 0; k 0]);

%% Coefficient Definition
c = [1 0 1 D 0 D]';
a = [0 0 1 0]';
f = [0 pres]';

%% Finite Element and Analytical Solutions
u = assempde(pdem,c,a,f);
numNodes = size(pdem.Mesh.Nodes,2);
figure
pdeplot(pdem, 'xydata', u(1:numNodes), 'contour', 'on');
title 'Transverse Deflection'

numNodes = size(pdem.Mesh.Nodes,2);
fprintf('Transverse deflection at plate center(PDE Toolbox) = %12.4e\n', min(u(1:numNodes,1)));
% compute analytical solution
wMax = -.0138*pres*len^4/(E*thick^3);
fprintf('Transverse deflection at plate center(analytical) = %12.4e\n', wMax);


