
h = 0.25;

gAm = 0.2;
mu = -2;  % position
sigma = 0.9;  % width
external_input = @(x, t) gAm * exp(-(x - mu).^2 / (2 * sigma^2));

psi = @(x) x.*exp(-x);  
x1 = 0.0;
x2 = x1 + fzero(@(x) x.*exp(-x) - h, 2.5); 

initial_profile = @(x) sign(x-x1).*psi(abs(x-x1)) + sign(x2-x).*psi(abs(x2-x));

%% Run Dynamics
[x, t, uHist] = finiteElementsCollocationMethod(h, external_input, initial_profile);

figure;
imagesc(t,x,uHist');
colorbar; 
title('Temperature plot');
xlabel('time coordinates');
ylabel('x coordinates');
axis xy;
