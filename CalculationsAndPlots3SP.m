% How to cite: Mohammed, M., Blasius, B., & Ryabov, A. (2021). 
% Coexistence patterns and diversity in a trait-based metacommunity 
% on an environmental gradient. bioRxiv.

%%
% Author: Mozzamil Mohammed
% ICBM, University of Oldenburg
% Last update of the code: December 2021
%%
% calcuations and plots of species biomass and diversity: 3 species
%%
clearvars;

set(0,'defaultAxesFontSize',20);
set(0, 'DefaultLineLineWidth', 3);

load('trade_off_curve_Results_of_3sp .mat')


Ds =  logspace(-3, 3, 100);  % different diffusion rates using log10 scale
ResRangeMin = linspace(20,0.5, 100 ); % minimal range

% local and total biomass of 3 species 

  tot_biomass = cell(length(ResRangeMin),length(Ds));   % total biomass
  local_biomass = cell(length(ResRangeMin),length(Ds)); % local biomass
  for i = 1:length(ResRangeMin)
      for j=1:length(Ds)
          TL = trade_off_curve_Results_of_3sp{i,j};          
          tot_biomass(i,j) = TL(1);
          local_biomass(i,j) = TL(2);
      end
  end
  %% 

 % total biomass for the competition outcomes of 3 species for 3 different
 % environmental variabilities
 
 competition_outcomes_3sp = cell(3,1);
 competition_outcomes_3sp{1,:} = tot_biomass(40,:);
 competition_outcomes_3sp{2,:} = tot_biomass(50,:);
 competition_outcomes_3sp{3,:} = tot_biomass(100,:);
 

 % getting the total biomass in a matrix
 
 M1 = NaN(length(Ds),k);
 M2 = NaN(length(Ds),k);
 M3 = NaN(length(Ds),k);
 for i=1:length(Ds)
     M1(i,:) = competition_outcomes_3sp{1}{i};
     M2(i,:) = competition_outcomes_3sp{2}{i};
     M3(i,:) = competition_outcomes_3sp{3}{i};
 end

 % calculating the regional diversity 
 
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

% local biomass for 3 competing species for 3 different environmental
% variability

 competition_outcomes_3sp_loc = cell(3,1);
 competition_outcomes_3sp_loc{1,:} =   local_biomass(40,:);
 competition_outcomes_3sp_loc{2,:} =   local_biomass(50,:);
 competition_outcomes_3sp_loc{3,:} =   local_biomass(100,:);
 

 % local biomass 
 
 M1_loc = cell(length(Ds),1);
 M2_loc = cell(length(Ds),1);
 M3_loc = cell(length(Ds),1);
 for i=1:length(Ds)
     M1_loc{i,1} = competition_outcomes_3sp_loc{1}{i};
     M2_loc{i,1} = competition_outcomes_3sp_loc{2}{i};
     M3_loc{i,1} = competition_outcomes_3sp_loc{3}{i};
 end
 

% calculating local diversity

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

figure(2);  

cm = summer;
cm = cm(end:-1:1, :);
colormap(cm);

M1(M1 < 1e-1) = NaN;
M2(M2 < 1e-1) = NaN;
M3(M3 < 1e-1) = NaN;


ax1 = subplot(2,3,1)
pcolor_central(Ds, [2 6 10], M1')  % ploting regional diversity
shading flat

ax1 = gca;
ax1.XScale = 'log';
xticks([])
yticks([2 6 10])

%xlabel('Dispersal rate')
ylabel('Trait value, R^*_{i,1}')
title('\Delta{S} = 16'); %\Delta{S} = 16.3636
text(ax1,0.0008,13.1,'(a)','FontSize',25)

% I CHANGED HERE
% new_ax = axes('position',get(ax1,'position'),'color','none');
% set(new_ax,'YAxisLocation','right', 'XAxisLocation','top')
% xticks(new_ax,[]);
% ylim([1 12])
% yticks([])
%ylabel(new_ax,'Trait value, R^*')


ax2 = subplot(2,3,2)
pcolor_central(Ds, [2 6 10], M2')  % ploting regional diversity
shading flat
ax2 = gca;
ax2.XScale = 'log';
xticks([])
yticks([])
%xlabel('Dispersal rate')
title(' \Delta{S} = 20'); %\Delta{S} = 20.3030
text(ax2,0.0008,13.1,'(b)','FontSize',25)

% I CHANGED HERE
% new_ax = axes('position',get(ax2,'position'),'color','none');
% set(new_ax,'YAxisLocation','right', 'XAxisLocation','top')
% xticks(new_ax,[]);

% ylim([1 12])
% yticks([])


ax3 = subplot(2,3,3)
pcolor_central(Ds, [2 6 10], M3')  % ploting regional diversity
shading flat
ax3 = gca;
ax3.XScale = 'log';
xticks([])
yticks([])
%xlabel('Dispersal rate')
h = colorbar;
ylabel(h,'Total biomass')
set(h,'Position',[0.91    0.6   0.0181    0.3057])
h.Ticks = [100 350 600]; %Create 8 ticks from zero to 1
h.TickLabels = [100 350 600];
title('\Delta{S} = 40');
text(ax3,0.0008,13.1,'(c)','FontSize',25)

%%%% I CHANGED HERE

% new_ax = axes('position',get(ax3,'position'),'color','none');
% set(new_ax,'YAxisLocation','right', 'XAxisLocation','top')
% xticks(new_ax,[]);
% ylim([1 12])
% yticks([linspace(2,10,k)])
% ylabel(new_ax,'Trait value, R^*_1')

ax4 = subplot(2,3,4)
plot(Ds,Regional_M1)
hold on
plot(Ds,local_diversity_M1_loc)
ax4 = gca;
ax4.XScale = 'log';
legend('Regional diversity','Local diversity','Location','northwest')
ylim([1 3]);
xticks([10^-3 10^0 10^3]);
yticks([1 2 3]);
xlabel('Dispersal rate');
ylabel('Diversity');
% title('\Delta{S} = 16'); %\Delta{S} = 16.3636
text(ax1,0.0008,-3.8,'(d)','FontSize',25)

 

ax5 = subplot(2,3,5)
plot(Ds,Regional_M2)
hold on
plot(Ds,local_diversity_M2_loc)
ax5 = gca;
ax5.XScale = 'log';
ylim([1 3]);
xticks([10^-3 10^0 10^3]);
yticks([]);
xlabel('Dispersal rate');
%title('(e) \Delta{S} = 20'); %\Delta{S} = 20.3030
text(ax2,0.0008,-3.8,'(e)','FontSize',25)

ax6 = subplot(2,3,6)
plot(Ds,Regional_M3)
hold on
plot(Ds,local_diversity_M3_loc)
ax6 = gca;
ax6.XScale = 'log';
ylim([1 3]);
xticks([10^-3 10^0 10^3]);
yticks([]);
xlabel('Dispersal rate');
% title('(f) \Delta{S} = 40');
text(ax3,0.0008,-3.8,'(f)','FontSize',25)


  
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



figure(1);
%subplot(3,2,3)
ax1 = subplot(2,2,1) % regional diversity 3 species
pcolor_central(Ds, linspace(1,40,Var_steps), RegionalH)  % ploting regional diversity
shading flat
%xlabel('Dispersal rate')
h = colorbar;
ylabel(h,'E_{reg}', 'rotation', 360)
%text(0.1,0.1,'E_{reg}')
set(h,'Position',[0.47    0.6   0.0181    0.3057])
h.Ticks = [1 1.5 2 2.5]; 
pos = get (h, 'position' );
h.Label.Position = [2.5  pos(1)+2.5]; % to change its position
% h = colorbar;
% ylabel(h, 'Regional diversity')
%xticks([])
%yticks([1 2 3])
yticks([1 10 20 30 40])
xticks([])
% c = colorbar;
% c.title = 'Regional Diversity';
ax1 = gca;
ax1.XScale = 'log';
%colorbar

ylabel('Resource variabiliy')
title('3 Species')
text(ax1,0.001,43.9, '(a)', 'FontSize', 25)

ax3 = subplot(2,2,3) % local diversity 3 species

pcolor_central(Ds, linspace(1,40,Var_steps), local_diversity) % ploting regional diversity 
shading flat
ylabel('Resource variability')
xlabel('Dispersal rate')
h = colorbar;
ylabel(h,'E_{loc}','rotation', 360)
set(h,'Position',[0.47    0.13   0.0181    0.3057])
h.Ticks = [1 1.5 2 2.5]; 
pos = get (h, 'position' );
h.Label.Position = [2.5  pos(1)+2.5]; % to change its position
% h = colorbar;
% ylabel(h, 'Local diversity')
yticks([1 10 20 30 40])
xticks([10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3])
% c = colorbar;
% c.title = 'Regional Diversity';
ax3 = gca;
ax3.XScale = 'log';
%caxis([0, 30])
%colorbar 
% title('3 Species')
text(ax3,0.001,43.9, '(c)', 'FontSize', 25)


