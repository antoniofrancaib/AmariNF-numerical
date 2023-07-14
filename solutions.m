h = 0.25;
constant = @(x,t) 0.0; 

psi = @(x) x.*exp(-x);
x1 = 0.0;
x2 = x1 + fzero(@(x) x.*exp(-x) - h, 1); 
initial_profile = @(x) sign(x-x1).*psi(abs(x-x1)) + sign(x2-x).*psi(abs(x2-x));

[x, t, uHist] = finiteElementsCollocationMethod(h, constant, initial_profile);

figure;
imagesc(t,x,uHist');
colorbar; 
title('Temperature plot');
xlabel('time coordinates');
ylabel('x coordinates');
axis xy;


h = 0.25;

Am = 1.0;  
k = pi/5; 
w = pi/15;  
phi = pi;  
sinusoidal = @(x, t) Am*sin(k*x + w*t + phi);

psi = @(x) x.*exp(-x);
x1 = 0.0;
x2 = x1 + fzero(@(x) x.*exp(-x) - h, 1); 
initial_profile = @(x) sign(x-x1).*psi(abs(x-x1)) + sign(x2-x).*psi(abs(x2-x));

[x, t, uHist] = finiteElementsCollocationMethod(h, sinusoidal, initial_profile);

figure;
imagesc(t,x,uHist');
colorbar; 
title('Temperature plot');
xlabel('time coordinates');
ylabel('x coordinates');
axis xy;
