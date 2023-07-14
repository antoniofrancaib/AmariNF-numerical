function [x, t, uHist] = finiteElementsCollocationMethod(h, external_input, initial_profile)

    %% Create the spatial grid
    nx = 2000; 
    L = pi; 
    x = linspace(-L,L,nx)'; 
    delta_x = (2.0*L)/(nx-1); 
    
    %% Synaptic kernel matrix 
    omega = @(x,y) (1.0 - abs(x-y)).*exp(-abs(x-y));
    rho = delta_x*[0.5; ones(nx-2,1); 0.5]; 

    W = zeros(nx,nx); 
    for i = 1:nx
        for j = 1:nx
            W(i,j) = omega(x(i),x(j))*rho(j);
        end 
    end 

    %% Firing Rate Function
    mu = 50; % condition: mu >> 1, sigmoid to Heaviside as mu --> inf; 
    sigmoid = @(u) 1.0./(1.0+exp(-mu*(u-h))); 

    Heaviside = @(u) (u-h >= 0); 
    %% RHS Function 
    N = @(t,u) -u + W*Heaviside(u) + external_input(x,t);

    %% Run Dynamics
    tspan = [0, 200];
    [t,uHist] = ode45(N,tspan,initial_profile(x));
end