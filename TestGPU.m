function [Dt]=TestGPU(N,Nsteps)
   
    t(1)=now;
    
    vvg = WaveGPU(N,Nsteps);
    t(2)=now;
    vvg = WaveCPU(N,Nsteps);
    t(3)=now;
    Dt(1)=t(2)-t(1);
    Dt(2)=t(3)-t(2);
end

function vvg = WaveCPU(N,Nsteps)
%% Solving 2nd Order Wave Equation Using Spectral Methods
% This example solves a 2nd order wave equation: utt = uxx + uyy, with u =
% 0 on the boundaries. It uses a 2nd order central finite difference in
% time and a Chebysehv spectral method in space (using FFT). 
%
% The code has been modified from an example in Spectral Methods in MATLAB
% by Trefethen, Lloyd N.

% Points in X and Y
x = cos(pi*(0:N)/N); % using Chebyshev points

% Send x to the GPU
% x = gpuArray(x);
y = x';

% Calculating time step
dt = 6/N^2;

% Setting up grid
[xx,yy] = meshgrid(x,y);

% Calculate initial values
vv = exp(-40*((xx-.4).^2 + yy.^2));
vvold = vv;

ii = 2:N;
index1 = 1i*[0:N-1 0 1-N:-1];
index2 = -[0:N 1-N:-1].^2;

% Sending data to the GPU
% dt = gpuArray(dt);
% index1 = gpuArray(index1);
% index2 = gpuArray(index2);

% Creating weights used for spectral differentiation 
W1T = repmat(index1,N-1,1);
W2T = repmat(index2,N-1,1);
W3T = repmat(index1.',1,N-1);
W4T = repmat(index2.',1,N-1);

WuxxT1 = repmat((1./(1-x(ii).^2)),N-1,1);
WuxxT2 = repmat(x(ii)./(1-x(ii).^2).^(3/2),N-1,1);
WuyyT1 = repmat(1./(1-y(ii).^2),1,N-1);
WuyyT2 = repmat(y(ii)./(1-y(ii).^2).^(3/2),1,N-1);

% Start time-stepping
n = 0;

while n < Nsteps
    [vv,vvold] = stepsolution(vv,vvold,ii,N,dt,W1T,W2T,W3T,W4T,...
                                 WuxxT1,WuxxT2,WuyyT1,WuyyT2);
    n = n + 1;
end


% Gather vvg back from GPU memory when done
vvg = (vv);

end

function vvg = WaveGPU(N,Nsteps)
%% Solving 2nd Order Wave Equation Using Spectral Methods
% This example solves a 2nd order wave equation: utt = uxx + uyy, with u =
% 0 on the boundaries. It uses a 2nd order central finite difference in
% time and a Chebysehv spectral method in space (using FFT). 
%
% The code has been modified from an example in Spectral Methods in MATLAB
% by Trefethen, Lloyd N.

% Points in X and Y
x = cos(pi*(0:N)/N); % using Chebyshev points

% Send x to the GPU
x = gpuArray(x);
y = x';

% Calculating time step
dt = 6/N^2;

% Setting up grid
[xx,yy] = meshgrid(x,y);

% Calculate initial values
vv = exp(-40*((xx-.4).^2 + yy.^2));
vvold = vv;

ii = 2:N;
index1 = 1i*[0:N-1 0 1-N:-1];
index2 = -[0:N 1-N:-1].^2;

% Sending data to the GPU
dt = gpuArray(dt);
index1 = gpuArray(index1);
index2 = gpuArray(index2);

% Creating weights used for spectral differentiation 
W1T = repmat(index1,N-1,1);
W2T = repmat(index2,N-1,1);
W3T = repmat(index1.',1,N-1);
W4T = repmat(index2.',1,N-1);

WuxxT1 = repmat((1./(1-x(ii).^2)),N-1,1);
WuxxT2 = repmat(x(ii)./(1-x(ii).^2).^(3/2),N-1,1);
WuyyT1 = repmat(1./(1-y(ii).^2),1,N-1);
WuyyT2 = repmat(y(ii)./(1-y(ii).^2).^(3/2),1,N-1);

% Start time-stepping
n = 0;

while n < Nsteps
    [vv,vvold] = stepsolution(vv,vvold,ii,N,dt,W1T,W2T,W3T,W4T,...
                                 WuxxT1,WuxxT2,WuyyT1,WuyyT2);
    n = n + 1;
end


% Gather vvg back from GPU memory when done
vvg = gather(vv);

end


function [vv,vvold] = stepsolution(vv,vvold,ii,N,dt,W1T,W2T,W3T,W4T,...
                                        WuxxT1,WuxxT2,WuyyT1,WuyyT2)
    V = [vv(ii,:) vv(ii,N:-1:2)];
    U = real(fft(V.')).';
    
    W1test = (U.*W1T).';
    W2test = (U.*W2T).';
    W1 = (real(ifft(W1test))).';
    W2 = (real(ifft(W2test))).';
    
    % Calculating 2nd derivative in x
    uxx(ii,ii) = W2(:,ii).* WuxxT1 - W1(:,ii).*WuxxT2;
    uxx([1,N+1],[1,N+1]) = 0;
    
    V = [vv(:,ii); vv((N:-1:2),ii)];
    U = real(fft(V));
    
    W1 = real(ifft(U.*W3T));
    W2 = real(ifft(U.*W4T));
    
    % Calculating 2nd derivative in y
    uyy(ii,ii) = W2(ii,:).* WuyyT1 - W1(ii,:).*WuyyT2;
    uyy([1,N+1],[1,N+1]) = 0;
    
    % Computing new value using 2nd order central finite difference in
    % time
    vvnew = 2*vv - vvold + dt*dt*(uxx+uyy);
    vvold = vv; vv = vvnew;

end