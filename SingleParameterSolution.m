% How to cite: Mohammed, M., Blasius, B., & Ryabov, A. (2021). 
% Coexistence patterns and diversity in a trait-based metacommunity 
% on an environmental gradient. bioRxiv.

%%
% Author: Mozzamil Mohammed
% ICBM, University of Oldenburg
% Last update of the code: December 2021

%%

% this function returns the solution of ModelEquations() after passing the
% parameter values, data

% to run this part, you must download and install the high performance computing 
% SANDIALS PACKAGE 
% (https://computing.llnl.gov/projects/sundials/sundials-software) or
% alternatively comment the CVODE part and uncomment the standard matlab
% ode solver, but then the compution would be very slow

function Results =  SingleParameterSolution(data)



% initial densities

%InitSpecies = linspace(1,20,data.k* data.n);
InitResource = linspace(1,20,2*data.n);
InitDensities = [data.InitSpecies InitResource]; % for species and 2 resources


% uncomment this part to use the standard Matlab ode solver

% initial and maximum time
t0 = 0;
%timeSpan = [0 200000];
%timeSpan = [0 10000];
% timeSpan = linspace(1,  10000, 3);
% % solutions
% options = odeset('RelTol',10^-3,'AbsTol',10^-3); 
% %options = odeset('RelTol',10^-3,'AbsTol',10^-3); 
% 
% [t,x] = ode23(@(t, x) ModelEquations(t,x, k,n,g_max, m, D, F, K1, K2, c1, c2, S1, S2),timeSpan,InitDensities(:),options);


% here is the numerical integration using CVODE
% Options for integration

options = CVodeSetOptions('UserData', data,...
                          'RelTol', 0,...
                          'AbsTol', 1.e-5,...
                          'LinearSolver','Dense');

% mondata.grph = false;
% options = CVodeSetOptions(options,...
%                           'MonitorFn',@CVodeMonitor,...
%                           'MonitorData',mondata);



                      
CVodeInit(@ModelEquations, 'BDF', 'Newton', t0, InitDensities(:), options);


[status, t,x] = CVode([1:10000], 'Normal');

% Free solver memory
CVodeFree;

% we slove until 10000 units of time and use the solution as an initial
% guess for finding the roots of the system using fsolve

fmodx= @(x) ModelEquations(t,x, data);


xfin = fsolve(fmodx, x(:,end));

   

x = [x xfin];
t = [t 1e6];

% reshaping the outputs
X = reshape(x', length(t), data.n, data.k +2);

% Display solver statistics
%Stats = CVodeGetStats


% getting the results 
Results = cell(2,1);
TotalBioMass = sum(X(end, :, 1:data.k), 2); % total biomass of species k in all patches
TotalBioMass = reshape(TotalBioMass(:), 1, data.k);
EachPatchEquilibrium = X(end, :, 1:data.k); % the equilibrium of species k in each patch n
EachPatchEquilibrium  = reshape(EachPatchEquilibrium(:), data.n, data.k);
Results{1} = TotalBioMass;
Results{2} = EachPatchEquilibrium;

end








