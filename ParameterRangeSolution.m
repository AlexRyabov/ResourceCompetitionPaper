% How to cite: Mohammed, M., Blasius, B., & Ryabov, A. (2021). 
% Coexistence patterns and diversity in a trait-based metacommunity 
% on an environmental gradient. bioRxiv.

%%
% Author: Mozzamil Mohammed
% ICBM, University of Oldenburg
% Last update of the code: December 2021

%%

% this function has no inputs and returns solutions for parameter ranges (dispersal rate and resource supply)
% for 3 and 15 species

function ParameterRangeSolution()
clearvars;

%% 3 species

Number_of_patches = 50;
Number_of_species = 3;
Var_steps = 100;
Diff_steps = 100;

k = Number_of_species;
n = Number_of_patches;


g_max = 1;      % maximum growth rate
m = 0.25;       % mortality rate
F = m;          % dilution rate

% defining two different trade-off curves
flag = 0;
if flag
    xi = linspace(-2, 2, Number_of_species); a0 = -2;
    R_star1 = 16./(1 + exp(-a0-xi)) + 1.8; % minimal requirement for resource 1
    R_star2 = 16./(1 + exp(-a0+xi)) + 1.8;  % minimal requirement for resource 2
    K1 =  R_star1*(g_max-m)/m;    % vector of half-saturation constants for species i on resource 1
    K2 =  R_star2*(g_max-m)/m;    % vector of half-saturation constants for species i on resource 2
    c1 =  0.05*R_star1;     % consumption vectors of species i on resource 1
    c2 =  0.05*R_star2;      % consumption vectors of species i on resource 2
else                          % vector of half-saturation constants for species i on resource 2
    R_star = linspace(2,10,Number_of_species); % minimal requirement for resource 1
    K1 =  R_star*(g_max-m)/m;    % vector of half-saturation constants for species i on resource 1
    K2 =  K1(end:-1:1);     % vector of half-saturation constants for species i on resource 2
    c1 =  0.05*R_star;     % consumption vectors of species i on resource 1
    c2 =  c1(end:-1:1);      % consumption vectors of species i on resource 2
end


% different diffusion rates using log10 scale
Ds =  logspace(-3, 3, Diff_steps);      


data.k = k;  
data.n = n;
data.g_max = g_max;
data.m = m;
%data.D = D;
data.F = F;
data.K1 = K1;
data.K2 = K2;
data.c1 = c1;
data.c2 = c2;
%data.S1 = S1;
%data.S2 = S2;

% Initial densities: uniform and random distribution of species along the
% spatial gradient

flag = 0;
if flag
    % random
    ind = randperm(data.k* data.n);
    InitSpecies_rand = linspace(1,20,data.k* data.n);
    data.InitSpecies = InitSpecies_rand(ind)
else 
    % uniform
    data.InitSpecies = linspace(1,20,data.k* data.n);
end

% spatial environmental variability 

ResRangeMin = linspace(20,0.5,Var_steps ); % minimal range
ResRangeMax = linspace(21.,40.5,Var_steps ); % maximal range

Res = NaN(Var_steps ,Number_of_patches);

FileName = ['trade_off_curve_Results_of_3sp .mat'] % name of the file where results to be saved
%Results_of_3sp
%trade_off_curve_Results_of_3sp
%trade_off_curve_Results_of_3sp

Calculate = 1;

if Calculate
    
    % random and uniform distribution of resources along the spatial
    % gradient
    trade_off_curve_Results_of_3sp = cell(length(ResRangeMin),length(Ds));
    Srandperm = [32 21 8 6 18 5 36 46 23 20 35 49 33 17 26 45 ...
        40 13 28 42 34 38 31 24 12 16 9 10 19 1 50 14 ...
        29 25 7 30 27 15 22 44 41 2 47 43 37 11 39 3 48 4];
    for i=1:length(ResRangeMin)
        ii = i
        Res(i,:) = linspace(ResRangeMin(i),ResRangeMax(i),Number_of_patches);
        FlagRandom = 0;
        if FlagRandom
            % random
            ind = Srandperm;
            S1 = Res(i,ind); % supply of resource 1
            S2 = S1(end:-1:1); % supply of resource 2
        else
            % uniform
            S1 = Res(i,:);      
            S2 = S1(end:-1:1);                
        end
        data.S1 = S1;
        data.S2 = S2;
        for j = 1:length(Ds)
            D = Ds(j);
            data.D = D;
            trade_off_curve_Results_of_3sp{i, j} = SingleParameterSolution(data);
        end
    end
    
    save(FileName, 'trade_off_curve_Results_of_3sp');
else
    load(FileName, '-mat');
end
 


%% 

% 15 species

Number_of_species = 15;
k = Number_of_species;


% defining two different trade-off curves
flag = 0;
if flag
    xi = linspace(-2, 2, Number_of_species); a0 = -2;
    R_star1 = 16./(1 + exp(-a0-xi)) + 1.8; % minimal requirement for resource 1
    R_star2 = 16./(1 + exp(-a0+xi)) + 1.8;  % minimal requirement for resource 2
    K1 =  R_star1*(g_max-m)/m;    % vector of half-saturation constants for species i on resource 1
    K2 =  R_star2*(g_max-m)/m;    % vector of half-saturation constants for species i on resource 2
    c1 =  0.05*R_star1;     % consumption vectors of species i on resource 1
    c2 =  0.05*R_star2;      % consumption vectors of species i on resource 2
else                          % vector of half-saturation constants for species i on resource 2
    R_star = linspace(2,10,Number_of_species); % minimal requirement for resource 1
    K1 =  R_star*(g_max-m)/m;    % vector of half-saturation constants for species i on resource 1
    K2 =  K1(end:-1:1);     % vector of half-saturation constants for species i on resource 2
    c1 =  0.05*R_star;     % consumption vectors of species i on resource 1
    c2 =  c1(end:-1:1);      % consumption vectors of species i on resource 2
end


data.k = k;  


% Initial densities: uniform and random distribution of species along the
% spatial gradient

flag = 0;
if flag
    % random
    ind = randperm(data.k* data.n);
    InitSpecies_rand = linspace(1,20,data.k* data.n);
    data.InitSpecies = InitSpecies_rand(ind)
else 
    % uniform
    data.InitSpecies = linspace(1,20,data.k* data.n);
end

% spatial environmental variability 

ResRangeMin = linspace(20,0.5,Var_steps ); % minimal range
ResRangeMax = linspace(21.,40.5,Var_steps ); % maximal range

Res = NaN(Var_steps ,Number_of_patches);

FileName = ['trade_off_curve_Results_of_15sp .mat']
%Results_of_15sp
%trade_off_curve_Results_of_15sp
%trade_off_curve_Results_of_15sp

Calculate = 1;

if Calculate
    % uniform and random distribution of resources along the gradient
    trade_off_curve_Results_of_15sp = cell(length(ResRangeMin),length(Ds));
    Srandperm = [32 21 8 6 18 5 36 46 23 20 35 49 33 17 26 45 ...
        40 13 28 42 34 38 31 24 12 16 9 10 19 1 50 14 ...
        29 25 7 30 27 15 22 44 41 2 47 43 37 11 39 3 48 4];
    for i=1:length(ResRangeMin)
        ii = i
        Res(i,:) = linspace(ResRangeMin(i),ResRangeMax(i),Number_of_patches);
        FlagRandom = 0;
        if FlagRandom
            % random
            ind = Srandperm;
            S1 = Res(i,ind);
            S2 = S1(end:-1:1);
        else
            % uniform
            S1 = Res(i,:);      % supply of resource 1
            S2 = S1(end:-1:1);                % supply of resource 2
        end                % supply of resource 2
        data.S1 = S1;
        data.S2 = S2;
        for j = 1:length(Ds)
            D = Ds(j);
            data.D = D;
            trade_off_curve_Results_of_15sp{i, j} = SingleParameterSolution(data);
        end
    end
    
    save(FileName, 'trade_off_curve_Results_of_15sp');
else
    load(FileName, '-mat');
end
%% 

end