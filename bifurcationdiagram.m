%% Cleaning
clear,  clc;

%% External Input Function
external_input = @(x, t) 0.0;

%% Bifurcation Analysis Validation 
h_min = 0.01;
h_max = 0.4;
nh = 30;

h_values = linspace(h_min, h_max, nh); 
delta_values = zeros(size(h_values)); 

for i = 1:length(h_values)
    
    h = h_values(i);

    % initial Profile Function for this h
    psi = @(x) x.*exp(-x);  
    x1 = 0.0;
    x2 = x1 + fzero(@(x) x.*exp(-x) - h, 2.5); 
    
    initial_profile = @(x) sign(x-x1).*psi(abs(x-x1)) + sign(x2-x).*psi(abs(x2-x));

    [x, t, uHist] = finiteElementsCollocationMethod(h, external_input, initial_profile);

    uFinal = uHist(end,:);

    % measure delta 
    [~, idx] = mink(abs(uFinal - h), 2);
    delta_values(i) = x(max(idx)) - x(min(idx));
end

%% Theoretical Bifurcation Diagram 
delta_min = 0;
delta_max = 8;
ndelta = 200;

psi = @(x) x.*exp(-x);
delta = linspace(delta_min, delta_max, ndelta); 

figure;
plot(psi(delta), delta);
hold on; 
plot(h_values, delta_values, 'o')
hold off; 
title('Bifurcation Diagram Validation');
xlabel('h');
ylabel('\Delta');