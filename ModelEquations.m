
% this function defines the resouce competition model described in section
% 2 of the paper
function [dxdt, flag, new_data] = ModelEquations(t,x, data)

k = data.k;  % number of species
n = data.n;  % number of patches
g_max = data.g_max; % maximal growth rate
m = data.m;  % mortality rate
D = data.D;  % diffusion rate
F = data.F;  % dilution rate
K1 = data.K1;% half-saturation constant for resource 1
K2 = data.K2;% half-saturation constant for resource 2
c1 = data.c1;% consumption rate for resource 1
c2 = data.c2; % consumption rate for resource 2
S1 = data.S1; % environmental supply of resource 1
S2 = data.S2; % environmental supply of resource 2

% pre-locations
dxdt = zeros(n, k + 2); % derivatives of species Ns and resources Rs

N = x(1:(k*n));  % N is the average density of the species k in patch n
N = reshape(N, n, k); % reshaping: a matrix has n rows and k columns 

R1 = x(k*n+1:((k+1)*n));  % resource 1


R2 = x((k+1)*n+1:end);    % resource 2



gr = g_max*min(R1(:)./(K1(1:k)+R1(:)),R2(:)./(K2(1:k)+R2(:))) .* N(:,:); % growth rate of all species in all patches





dxdt(1,1:k) = gr(1,:) - m*N(1,:) + D*(N(2,:) - N(1,:)); % overall growth rate of species on the left boundary and can only disperse to the right cell



dxdt(2:n-1,1:k) = gr(2:n-1,:) - m*N(2:n-1,:) + D*(N(1:n-2,:) - 2*N(2:n-1,:) + N(3:n,:)); % overall growth rate of species i (i =2:n-1) and can disperse to its left and right cells

dxdt(n,1:k) = gr(n,:) - m*N(n,:) + D*(N(n-1,:) - N(n,:));          % overall growth rate of species n on the right boundary and can only disperse to the left cell

dxdt(:,end-1) = F*(S1(:)-R1(:)) - (c1 * gr')'; % overall growth rate of resource 1
dxdt(:,end)   = F*(S2(:)-R2(:)) - (c2 * gr')'; % overall growth rate of resource 2

flag = 0;
new_data = [];

dxdt = dxdt(:);
        



end