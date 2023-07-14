%% Finite-Elements Collocation Method

%% Cleaning
clear,  clc;

%% Create the spatial grid
nx = 2000; 
L = pi; % coordinates bounds of the bar 
x = linspace(-L,L,nx)'; 
delta_x = (2.0*L)/(nx-1); % space differential 

%% Synaptic connectivity functions
omega = @(x,y) (1.0 - abs(x-y)).*exp(-abs(x-y));

%% Discretization: Synaptic kernel matrix 
rho = delta_x*[0.5; ones(nx-2,1); 0.5]; 

W = zeros(nx,nx); 
for i = 1:nx
    for j = 1:nx
        W(i,j) = omega(x(i),x(j))*rho(j);
    end 
end 

% observation: isequal(W(2:end-1, 2:end-1), W(2:end-1, 2:end-1)') - true

%% Firing Rate: use Sigmoid or Heaviside 
h = 0.25; % firing threshold; bifurcation occurs at h = 1/exp(1)
mu = 50; % condition: v >> 1, sigmoid to H as v --> inf; 
sigmoid = @(u) 1.0./(1.0+exp(-mu*(u-h))); 
H = @(u) (u-h >= 0); % heaviside function

%% External Input Function 
constant = @(x,t) 0.0; 

% Pulse parameters
t1 = 30; 
t2 = 100; 

% Sinusoidal function for traveling waves 
Am = 1.0;  % amplitude
k = 2*pi/10; % spatial frequency
w = pi/10;  % temporal frequency
phi = pi;  % phase
sinusoidal = @(x, t) Am*sin(k*x + w*t + phi);

% Use this for: replicate flies experiment
gAm = 0.2;
mu = -2;  % position
sigma = 1.5;  % width
pulse = @(x, t) gAm * exp(-(x - mu).^2 / (2 * sigma^2));%.* (t >= t1 & t <= t2);

%% Initial profile functions
A = 0.0; % height of the bump
alpha = 5.0; % 1/alpha proportional to width of the bump
phi = 1;
initial_profile = @(x) A ./ (cosh(alpha*(x - phi))).^2;

%% RHS Function 
N = @(t,u) -u + W*sigmoid(u) + constant(x,t);

% Use this for: Bifurcation Analysis Theory Validation 
psi = @(x) x.*exp(-x);
x1 = 0.0;
x2 = x1 + fzero(@(x) x.*exp(-x) - h, 1); 

theoretical_stable_bump = sign(x-x1).*psi(abs(x-x1)) + sign(x2-x).*psi(abs(x2-x));

%% Time step and solve the system of ODEs 
tspan = [0, 100];
[t,uHist] = ode45(N,tspan,theoretical_stable_bump);

%% Plot the approximation 
figure;
mesh(t,x,uHist')
title('Approximation to the Amari neural field model')
xlabel('time coordinates')
ylabel('x coordinates')
zlabel('Neuron Voltage')

figure;
imagesc(t,x,uHist');
colorbar; 
title('Temperature plot');
xlabel('time coordinates');
ylabel('x coordinates');
axis xy;


