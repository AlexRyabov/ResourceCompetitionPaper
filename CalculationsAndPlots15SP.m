% How to cite: Mohammed, M., Blasius, B., & Ryabov, A. (2021). 
% Coexistence patterns and diversity in a trait-based metacommunity 
% on an environmental gradient. bioRxiv.

%%
% Author: Mozzamil Mohammed
% ICBM, University of Oldenburg
% Last update of the code: December 2021
%%
% calcuations and plots of species biomass and diversity: 15 species
%%
clearvars;

set(0,'defaultAxesFontSize',20);
set(0, 'DefaultLineLineWidth', 3);

load('trade_off_curve_Results_of_15sp .mat')

Ds =  logspace(-3, 3, 100);  % different diffusion rates using log10 scale
ResRangeMin = linspace(20,0.5, 100 ); % minimal range

% local and total biomass of 15 species 

  tot_biomass = cell(length(ResRangeMin),length(Ds));   % total biomass
  local_biomass = cell(length(ResRangeMin),length(Ds)); % local biomass
  for i = 1:length(ResRangeMin)
      for j=1:length(Ds)
          TL = trade_off_curve_Results_of_15sp{i,j};          
          tot_biomass(i,j) = TL(1);
          local_biomass(i,j) = TL(2);
      end
  end

 % total biomass of 15 species for 3 different environmental variabilities
 
 competition_outcomes_3sp = cell(3,1);
 competition_outcomes_3sp{1,:} = tot_biomass(20,:);
 competition_outcomes_3sp{2,:} = tot_biomass(60,:);
 competition_outcomes_3sp{3,:} = tot_biomass(100,:);
 
% getting the the biomass in a matrix form

 M1 = NaN(length(Ds),k);
 M2 = NaN(length(Ds),k);
 M3 = NaN(length(Ds),k);
 for i=1:length(Ds)
     M1(i,:) = competition_outcomes_3sp{1}{i};
     M2(i,:) = competition_outcomes_3sp{2}{i};
     M3(i,:) = competition_outcomes_3sp{3}{i};
 end

% calculating regional diversity

Regional_M1 = [];  % reginal diversity
Regional_M2 = [];  % reginal diversity
Regional_M3 = [];  % reginal diversity

for i=1:length(Ds)
    M1_D = M1(i,:);
    M2_D = M2(i,:);
    M3_D = M3(i,:);
    Sum_Of_M1 = sum(M1_D);  % sum of total biomass in all patches for each D
    Rel_Abundance1 = M1_D./Sum_Of_M1;              % relative abundance of each species in all patches
    Sum_Of_M2 = sum(M2_D);  % sum of total biomass in all patches for each D
    Rel_Abundance2 = M2_D./Sum_Of_M2;              % relative abundance of each species in all patches
    Sum_Of_M3 = sum(M3_D);  % sum of total biomass in all patches for each D
    Rel_Abundance3 = M3_D./Sum_Of_M3;              % relative abundance of each species in all patches
    %LoG(i,:) = log(Rel_Abundance(i,:)+ 0.1e-15);            % ln(Relative abundance)
    ind1 = Rel_Abundance1 > 0;
    if max(ind1) > 0
        Regional_M1(end+1)  = exp(-sum( Rel_Abundance1(ind1) .* log(Rel_Abundance1(ind1)))); % Shannon index
    else
        Regional_M1 = 0;
    end
    
    ind2 = Rel_Abundance2 > 0;
    if max(ind2) > 0
        Regional_M2(end+1)  = exp(-sum( Rel_Abundance2(ind2) .* log(Rel_Abundance2(ind2)))); % Shannon index
    else
        Regional_M2 = 0;
    end
    
    ind3 = Rel_Abundance3 > 0;
    if max(ind3) > 0
        Regional_M3(end+1)  = exp(-sum( Rel_Abundance3(ind3) .* log(Rel_Abundance3(ind3)))); % Shannon index
    else
        Regional_M3 = 0;
    end
    
end
%% 




%% 

%% 

% local biomass of 15 species for 3 different environmanral variabilities

 competition_outcomes_3sp_loc = cell(3,1);
 competition_outcomes_3sp_loc{1,:} =   local_biomass(20,:);
 competition_outcomes_3sp_loc{2,:} =   local_biomass(60,:);
 competition_outcomes_3sp_loc{3,:} =   local_biomass(100,:);
 
 

 M1_loc = cell(length(Ds),1);
 M2_loc = cell(length(Ds),1);
 M3_loc = cell(length(Ds),1);
 for i=1:length(Ds)
     M1_loc{i,1} = competition_outcomes_3sp_loc{1}{i};
     M2_loc{i,1} = competition_outcomes_3sp_loc{2}{i};
     M3_loc{i,1} = competition_outcomes_3sp_loc{3}{i};
 end
 
 % species while changing evironmental variability for 4 diffusion rates
 
 species_var_total = cell(1,4);
 species_var_total{1} =   tot_biomass(:,1);
 species_var_total{2} =   tot_biomass(:,20); %40
 species_var_total{3} =   tot_biomass(:,50); %80
 species_var_total{4} =   tot_biomass(:,100);
 
 M1_total_Var = NaN(length(ResRangeMin),k);
 M2_total_Var = NaN(length(ResRangeMin),k);
 M3_total_Var = NaN(length(ResRangeMin),k);
 M4_total_Var = NaN(length(ResRangeMin),k);
 for i=1:length(ResRangeMin)
     M1_total_Var(i,:) =  species_var_total{1}{i};
     M2_total_Var(i,:) =  species_var_total{2}{i};
     M3_total_Var(i,:) =  species_var_total{3}{i};
     M4_total_Var(i,:) =  species_var_total{4}{i};
 end
 
 M1_total_Var( M1_total_Var < 1e-1) = NaN;
 M2_total_Var( M2_total_Var < 1e-1) = NaN;
 M3_total_Var( M3_total_Var < 1e-1) = NaN;
 M4_total_Var( M4_total_Var < 1e-1) = NaN;
 
 % species distributions in patches for a fixed variability (default 40 = M3_loc) and for 4
 % diffusion rates: I HAVE CHANGED M3 TO M1 TO GET SUPPLEMENTRY FIGURE
 
 species_distribution_D1 =  M3_loc{1};
 species_distribution_D2 =  M3_loc{20};
 species_distribution_D3 =  M3_loc{50};
 species_distribution_D4 =  M3_loc{100}; 
 
 species_distribution_D1(species_distribution_D1 < 1e-1) = NaN;
 species_distribution_D2(species_distribution_D2 < 1e-1) = NaN;
 species_distribution_D3(species_distribution_D3 < 1e-1) = NaN;
 species_distribution_D4(species_distribution_D4 < 1e-1) = NaN;
 
% local diversity  
local_diversity_M1_loc = [];  % local diversity
local_diversity_M2_loc = [];  % local diversity
local_diversity_M3_loc = [];  % local diversity
for i=1:length(Ds)
    loc1 =  M1_loc{i};
    loc1(loc1<0) = 0;
    Sum_Of_loc1 = sum(loc1, 2);
    Rel_Abundance_loc1 = loc1./ repmat(Sum_Of_loc1, 1, k);
    
    local_diversity_M1_loc(end+1) = mean(exp(-sum((Rel_Abundance_loc1) .* log(Rel_Abundance_loc1 + 1e-18), 2)),1);
    
    loc2 =  M2_loc{i};
    loc2(loc2<0) = 0;
    Sum_Of_loc2 = sum(loc2, 2);
    Rel_Abundance_loc2 = loc2./ repmat(Sum_Of_loc2, 1, k);
    
    local_diversity_M2_loc(end+1) = mean(exp(-sum((Rel_Abundance_loc2) .* log(Rel_Abundance_loc2 + 1e-18), 2)),1);
    
    loc3 =  M3_loc{i};
    loc3(loc3<0) = 0;
    Sum_Of_loc3 = sum(loc3, 2);
    Rel_Abundance_loc3 = loc3./ repmat(Sum_Of_loc3, 1, k);
    
    local_diversity_M3_loc(end+1) = mean(exp(-sum((Rel_Abundance_loc3) .* log(Rel_Abundance_loc3 + 1e-18), 2)),1);
end

%% 

%% 

figure(3); 
cm = summer;
cm = cm(end:-1:1, :);
colormap(cm);

M1(M1 < 1e-1) = NaN;
M2(M2 < 1e-1) = NaN;
M3(M3 < 1e-1) = NaN;


ax1 = subplot(2,3,1)
pcolor_central(Ds, linspace(2,10,k), M1')  % ploting regional diversity
shading flat
ax1 = gca;
ax1.XScale = 'log';
xticks([])
yticks([2 4 6 8 10])
%xlabel('Dispersal rate')
ylabel('Trait value, R^*_{i,1}')
title('\Delta{S} = 8'); %\Delta{S} = 8.4848
text(ax1,0.0008,17,'(a)','FontSize',25)

% I CHANGED HERE
% new_ax = axes('position',get(ax1,'position'),'color','none');
% set(new_ax,'YAxisLocation','right', 'XAxisLocation','top')
% xticks(new_ax,[]);
% ylim(new_ax,[2 10])
% yticks(new_ax,[])
  
ax2 = subplot(2,3,2)
pcolor_central(Ds, 1:k, M2')  % ploting regional diversity
shading flat
ax2 = gca;
ax2.XScale = 'log';
xticks([])
yticks([])
%xlabel('Dispersal rate')
title('\Delta{S} =  24'); %\Delta{S} =  24.2424
text(ax2,0.0008,17,'(b)','FontSize',25)


%  I CHANGED HERE
% new_ax = axes('position',get(ax2,'position'),'color','none');
% set(new_ax,'YAxisLocation','right', 'XAxisLocation','top')
% xticks(new_ax,[]);
% ylim(new_ax,[2 10])
% yticks(new_ax,[])
%ylabel(new_ax,'Trait value, R^*')

ax3 = subplot(2,3,3)
pcolor_central(Ds, 1:k, M3')  % ploting regional diversity
shading flat
ax3 = gca;
ax3.XScale = 'log';
xticks([])
yticks([])
%xlabel('Dispersal rate')
%h = colorbar;
%ylabel(h,'Total biomass')
h = colorbar;
ylabel(h,'Total biomass')
set(h,'Position',[0.91    0.6   0.0181    0.3057])
h.Ticks = [100 350 600];  
h.TickLabels = [100 350 600];

title('\Delta{S} = 40');
text(ax3,0.0008,17,'(c)','FontSize',25)
% 
% I CHANGED HERE
% new_ax = axes('position',get(ax3,'position'),'color','none');
% set(new_ax,'YAxisLocation','right', 'XAxisLocation','top')
% xticks(new_ax,[]);
% ylim(new_ax,[2 10])
% yticks(new_ax,[2    4    6  8   10])
% ylabel(new_ax,'Trait value, R^*_1')


ax4 =subplot(2,3,4)
plot(Ds,Regional_M1)
hold on
plot(Ds,local_diversity_M1_loc)
ax4 = gca;
ax4.XScale = 'log';
legend('Regional diversity','Local diversity','Location','northwest')
ylim([1 15]);
xticks([10^-3 10^0 10^3]);
yticks([1 5 10 15])
xlabel('Dispersal rate');
ylabel('Diversity');
% title('\Delta{S} = 8'); %\Delta{S} = 8.4848
text(ax1,0.0008,-4.3,'(d)','FontSize',25)


ax5 = subplot(2,3,5)
plot(Ds,Regional_M2)
hold on
plot(Ds,local_diversity_M2_loc)
ax5 = gca;
ax5.XScale = 'log';
ylim([1 15]);
xticks([10^-3 10^0 10^3]);
yticks([]);
xlabel('Dispersal rate');
%title('(e) \Delta{S} =  24'); %\Delta{S} =  24.2424
text(ax2,0.0008,-4.3,'(e)','FontSize',25)

subplot(2,3,6)
plot(Ds,Regional_M3)
hold on
plot(Ds,local_diversity_M3_loc)
ax = gca;
ax.XScale = 'log';
ylim([1 15]);
xticks([10^-3 10^0 10^3]);
yticks([]);
xlabel('Dispersal rate');
% title('(f) \Delta{S} = 40');
text(ax3,0.0008,-4.3,'(f)','FontSize',25)





figure(4);

cm = summer;
cm = cm(end:-1:1, :);
colormap(cm);


ax1 = subplot(2,2,1) 

pcolor_central(1:n, linspace(2,10,k),  species_distribution_D1')   
shading flat
xlim([1 n])
xticks([])
ylim([2 10])
yticks([2 4 6 8 10])
% 
% h = colorbar;
% ylabel(h, 'Regional diversity')
title('d = 0.001')
text(ax1, 0.5,10.8, '(a)', 'FontSize', 25)
%xlabel('Patches');
ylabel('Trait value, R^*_{i,1}');
% 
% new_ax = axes('position',get(ax1,'position'),'color','none');
% set(new_ax,'YAxisLocation','right', 'XAxisLocation','top')
% xticks(new_ax,[]);
% ylim(new_ax,[2 10])
% yticks(new_ax,[])
%ylabel(new_ax,'Trait value, R^*')


ax2 = subplot(2,2,2) 
pcolor_central(1:n, linspace(2,10,k),  species_distribution_D2')   
shading flat
xlim([1 n])
ylim([2 10])
yticks([2 4 6 8 10])
xticks([])
yticks([])
title('d = 0.01') %D = 0.0142
text(ax2, 0.5,10.8, '(b)', 'FontSize', 25)

h = colorbar;
ylabel(h,'Biomass')
set(h,'Position',[0.91    0.13   0.0181    0.3057])
h.Ticks = [10 20 30 40];  
% h.TickLabels = [100 350 600];
% 
% new_ax = axes('position',get(ax2,'position'),'color','none');
% set(new_ax,'YAxisLocation','right', 'XAxisLocation','top')
% xticks(new_ax,[]);
% ylim(new_ax,[2 10])
% yticks(new_ax,[2    4    6  8   10])
% ylabel(new_ax,'Trait value, R^*_1')


ax3 = subplot(2,2,3) 
pcolor_central(1:n, linspace(2,10,k),  species_distribution_D3')   
shading flat
xlim([1 n])
ylim([2 10])
yticks([2 4 6 8 10])
% h = colorbar;
% ylabel(h, 'Regional diversity')
xticks([1  10  20  30  40  50])
title('d = 1') %D = 0.9326
text(ax3, 0.5,10.8, '(c)', 'FontSize', 25)

xlabel('Patch number');
ylabel('Trait value, R^*_{i,1}');
% 
% new_ax = axes('position',get(ax3,'position'),'color','none');
% set(new_ax,'YAxisLocation','right', 'XAxisLocation','top')
% xticks(new_ax,[]);
% ylim(new_ax,[2 10])
% yticks(new_ax,[])
%ylabel(new_ax,'Trait value, R^*')

ax4 = subplot(2,2,4) 
pcolor_central(1:n, linspace(2,10,k),  species_distribution_D4')   
shading flat
xlim([1 n])
ylim([2 10])
yticks([2 4 6 8 10])
% h = colorbar;
% ylabel(h, 'Regional diversity')
xticks([1  10  20  30  40  50])
yticks([])
title('d = 1000')
text(ax4, 0.5,10.8, '(d)', 'FontSize', 25)

xlabel('Patch number');
% 
% new_ax = axes('position',get(ax4,'position'),'color','none');
% set(new_ax,'YAxisLocation','right', 'XAxisLocation','top')
% xticks(new_ax,[]);
% ylim(new_ax,[2 10])
% yticks(new_ax,[2    4    6  8   10])
% ylabel(new_ax,'Trait value, R^*_1') 



figure(5);

cm = summer;
cm = cm(end:-1:1, :);
colormap(cm);


ax1 = subplot(2,2,1) 
pcolor_central(linspace(1,40,100), linspace(2,10,k),  M1_total_Var')   
shading flat
xlim([1 40])
ylim([2 10])
yticks([2 4 6 8 10])
ylabel('Trait value, R^*_{i,1}');

% h = colorbar;
% ylabel(h, 'Regional diversity')
xticks([])

title('d = 0.001')
text(ax1, 0.5,10.8, '(a)', 'FontSize', 25)

% 
% new_ax = axes('position',get(ax1,'position'),'color','none');
% set(new_ax,'YAxisLocation','right', 'XAxisLocation','top')
% xticks(new_ax,[]);
% ylim(new_ax,[2 10])
% yticks(new_ax,[])
%ylabel(new_ax,'Trait value, R^*')

ax2 = subplot(2,2,2) 
pcolor_central(linspace(1,40,100), linspace(2,10,k),  M2_total_Var')   
shading flat
xlim([1 40])
ylim([2 10])
yticks([2 4 6 8 10])

xticks([])
yticks([])
title('d = 0.01') %D = 0.231
text(ax2, 0.5,10.8, '(b)', 'FontSize', 25)


h = colorbar;
ylabel(h,'Total biomass')
set(h,'Position',[0.91    0.13   0.0181    0.3057])
h.Ticks = [500 1000 1500 2000]; 

% new_ax = axes('position',get(ax2,'position'),'color','none');
% set(new_ax,'YAxisLocation','right', 'XAxisLocation','top')
% xticks(new_ax,[]);
% ylim(new_ax,[2 10])
% yticks(new_ax,[2    4    6  8   10])
% ylabel(new_ax,'Trait value, R^*_1')

ax3 = subplot(2,2,3) 
pcolor_central(linspace(1,40,100), linspace(2,10,k),  M3_total_Var')   
shading flat
xlim([1 40])
ylim([2 10])
yticks([2 4 6 8 10])
% h = colorbar;
% ylabel(h, 'Regional diversity')
xticks([1  10  20  30  40])
xlabel('Resource variability');
ylabel('Trait value, R^*_{i,1}');

title('d = 1') %D = 61.3591
text(ax3, 0.5,10.8, '(c)', 'FontSize', 25)



% 
% new_ax = axes('position',get(ax3,'position'),'color','none');
% set(new_ax,'YAxisLocation','right', 'XAxisLocation','top')
% xticks(new_ax,[]);
% ylim(new_ax,[2 10])
% yticks(new_ax,[])
%ylabel(new_ax,'Trait value, R^*')

ax4 = subplot(2,2,4) 
pcolor_central(linspace(1,40,100), linspace(2,10,k),  M4_total_Var')   
shading flat
xlim([1 40])
ylim([2 10])
yticks([2 4 6 8 10])
% h = colorbar;
% ylabel(h, 'Regional diversity')
xticks([1  10  20  30  40])
yticks([])
xlabel('Resource variability');
title('d = 1000')
text(ax4, 0.5,10.8, '(d)', 'FontSize', 25)



% 
% new_ax = axes('position',get(ax4,'position'),'color','none');
% set(new_ax,'YAxisLocation','right', 'XAxisLocation','top')
% xticks(new_ax,[]);
% ylim(new_ax,[2 10])
% yticks(new_ax,[2    4    6  8   10])
% ylabel(new_ax,'Trait value, R^*_1')




  
RegionalH = NaN(length(ResRangeMin),length(Ds));  % reginal diversity
for i=1:length(ResRangeMin)
    for j=1:length(Ds)
        Mi =  tot_biomass{i,j};
        Mi(Mi<0)=0;
        Sum_Of_Mi = sum(Mi);  % sum of total biomass in all patches for each D
        Rel_Abundance = Mi./Sum_Of_Mi;              % relative abundance of each species in all patches
        %LoG(i,:) = log(Rel_Abundance(i,:)+ 0.1e-15);            % ln(Relative abundance)
        ind = Rel_Abundance > 0;
        if max(ind) > 0
            RegionalH(i,j)  = exp(-sum( Rel_Abundance(ind) .* log(Rel_Abundance(ind)))); % Shannon index
        else
            RegionalH(i,j) = 0;
        end
        
        %      RegionalH(i)  = exp(-sum(mean(Rel_Abundance,1) .* log(mean(Rel_Abundance,1) + 0.1e-18))); % Shannon index
        %      if RegionalH(i) < 1e-6
        %          RegionalH(i) = 0;
        %      end
    end
end
 


local_diversity = NaN(length(ResRangeMin),length(Ds));  % local diversity
for i=1:length(ResRangeMin)
    for j=1:length(Ds)
        Ej =    local_biomass{i,j};
        Ej(Ej<0) = 0;
        Sum_Of_E = sum(Ej, 2);
        Rel_Abundance = Ej./ repmat(Sum_Of_E, 1, k);
        
        local_diversity(i,j) = mean(exp(-sum((Rel_Abundance) .* log(Rel_Abundance + 1e-18), 2)),1);
        
    end
end

%  LocalDiversity = mean(ESN,1); 
% 
figure(7);
ax1 = subplot(2,2,1);
plot(ax1, linspace(1,40,Var_steps), RegionalH(:,1), 'b')
xlim([1 40]);
ylim([1 15]);
yticks(ax1, [1 5 10 15]);
xticks(ax1, []);
ylabel('Diversity');
title('d = 0.001')
text(ax1, 0.5,16.5, '(a)', 'FontSize', 25)
hold on 
plot(ax1, linspace(1,40,Var_steps), local_diversity(:,1), 'r')
hold off

ax2 = subplot(2,2,2);
plot(ax2, linspace(1,40,Var_steps), RegionalH(:,20), 'b')
xlim([1 40]);
ylim([1 15]);
yticks(ax2, []);
xticks(ax2, []);
title('d = 0.01') %D = 0.231
text(ax2, 0.5,16.5, '(b)', 'FontSize', 25)
hold on 
plot(ax2, linspace(1,40,Var_steps), local_diversity(:,20), 'r')
hold off

ax3 = subplot(2,2,3);
plot(ax3, linspace(1,40,Var_steps), RegionalH(:,50), 'b')
xlim([1 40]);
ylim([1 15]);
yticks(ax3, [1 5 10 15]);
xticks(ax3, [1 10 20 30 40]);
xlabel('Resource variability');
ylabel('Diversity');
title('d = 1') %D = 61.3591
text(ax3, 0.5,16.5, '(c)', 'FontSize', 25)
hold on 
plot(ax3, linspace(1,40,Var_steps), local_diversity(:,50), 'r')
hold off

ax4 = subplot(2,2,4);
plot(ax4, linspace(1,40,Var_steps), RegionalH(:,100), 'b')
xlim([1 40]);
ylim([1 15]);
yticks(ax4, []);
xticks(ax4, [1 10 20 30 40]);
xlabel('Resource variability');
title('d = 1000')
text(ax4, 0.5,16.5, '(d)', 'FontSize', 25)

hold on 
plot(ax4, linspace(1,40,Var_steps), local_diversity(:,100), 'r')
hold off



figure(1);
%subplot(3,2,3)
ax2 = subplot(2,2,2) % regional diversity 15 species
pcolor_central(Ds, linspace(1,40,Var_steps), RegionalH)  % ploting regional diversity
shading flat
%xlabel('Dispersal rate')

h = colorbar;
ylabel(h,'E_{reg}', 'rotation', 360)
set(h,'Position',[0.91    0.6   0.0181    0.3057])
h.Ticks = [1 4 7 10 12]; 
pos = get (h, 'position' );
h.Label.Position = [2.4  pos(1)+14]; % to change its position
%xticks([])
%yticks([1 2 3])
yticks([])
xticks([])
% c = colorbar;
% c.title = 'Regional Diversity';
ax2 = gca;
ax2.XScale = 'log';
%colorbar

%ylabel('Resource supply variabiliy')
title('15 Species')
text(ax2,0.001,43.9, '(b)', 'FontSize', 25)

ax4 = subplot(2,2,4) % local diversity 15 species

pcolor_central(Ds, linspace(1,40,Var_steps), local_diversity) % ploting regional diversity 
shading flat
%ylabel('Resource variability')
xlabel('Dispersal rate')

h = colorbar;
ylabel(h,'E_{loc}','rotation',360)
set(h,'Position',[0.91    0.13   0.0181    0.3057])
h.Ticks = [1 3 5]; 
pos = get (h, 'position' );
h.Label.Position = [2.4  pos(1)+5.5]; % to change its position
yticks([])
xticks([10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3])
% c = colorbar;
% c.title = 'Regional Diversity';
ax4 = gca;
ax4.XScale = 'log';
%caxis([0, 30])
%colorbar 
% title('15 Species')
text(ax4,0.001,43.9, '(d)', 'FontSize', 25)
