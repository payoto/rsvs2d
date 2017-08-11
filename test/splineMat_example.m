%% splineMat example

order = 3;
Npts = 300;  % Number of points on aerofoil
Ncp = 10;    % Number of spline control points
[x,y] = naca(0,0,12,Npts);

% Input curve knots, u
%  (Parameterised by arc distance)
du = sqrt( diff(x).^2 + diff(y).^2 );
u = [0, cumsum(du)]/sum(du);


% Control point knots, v
%  (Choose arbitrary cosine distribution)
v1 = linspace(0,0.5,1+(Ncp-order)/2);
v1 = 0.5*sin(0:(pi/(1+Ncp-order)):pi/2);
%v = [v1(1:end-1), 1-flip(v1)];                     % No  fixed LE
v = [v1(1:end-1), ones(1,2)*v1(end) , 1-flip(v1)];  % Use fixed LE
v(end) = 1;


% Generate spline matrix
%  (Evaluates spline defined by knots v at knots u)
Nmat = splineMat(u, v, order);


% Remove repeated TE point in control points
Nmat(:,1) = Nmat(:,1) + Nmat(:,end);
Nmat = Nmat(:,1:end-1);

% Evaluation of last point equal to evaluation of first point (TE point)
Nmat(end,:) = Nmat(1,:);


% Inverse spline fitting
%  (Calculate control points from curve)
CP = Nmat\[x', y'];

figure;hold on;grid on;
plot(x,y);
plot(CP(:,1),CP(:,2),'o-');
