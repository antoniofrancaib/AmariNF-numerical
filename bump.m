h = 0.2;

A = 0.23; % height of the bump
alpha = 5.0; % 1/alpha proportional to width of the bump
initial_profile = A ./ (cosh(alpha*x)).^2;

psi = @(x) x.*exp(-x);
x1 = 0.0;
x2 = fzero(@(x) x.*exp(-x) - h, 0); 

theoretical_stable_bump = sign(x-x1).*psi(abs(x-x1)) + sign(x2-x).*psi(abs(x2-x));

plot(x,theoretical_stable_bump)
hold on; 
plot(x,initial_profile)
% yline(h)
% xline(x1); 
% xline(x2);
hold off; 
legend('Theoretical Stable Bump','Arbitrary Initial Profile', 'y=h','x=x1','x=x2');


% 2nd Option:
mu = -1;  % position
sigma = 0.2;  % width
bump_profile = exp(-(x - mu).^2 / (2 * sigma^2)); 