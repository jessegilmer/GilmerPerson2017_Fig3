function [Experiment_bin] = GCL_Main()
%% Created 11/23/16 by Gilmer, J.
% This is the main script. Control of the subroutines should start from
% here. Alternately, you can "load" in a system and run this script from
% the experimental section, if you don't want to make a new system every
% time.
%
% Here are the main variables:
%      GCL_Space - This is the 3D Size of the system the model exists in, the units are microns.
%                  It shold be a 3 variable vector, i.e. [1E3 1E3 1E3] makes a millimeter cube.
%
%      MFr_Density - This determines density of Mossy Fiber Rosettes (MFr) per cubic millimeter.
%                  The default value is 98,800 per Palkovits, 1972.
%
%      MFr_radius- This is the radius of an average MFr
%                  5 um is my conservative estimate, from several sources
%                  which vary around this point. See Sultan, 2001.
%
%      GC_Density -  This determines density of Granule Cells (GC) per cubic millimeter.
%                    The default value is  2,600,000 per Palkovits, 1971.
%
%      GC_MFr_Convergence - This is the average convergence onto a MFr by GCs.
%                  The default value is 56, back calculated from the 1:14
%                  population ratio.
%
%      GC_MFr_Divergence- This is the average divergence onto MFrs by a GC.
%                  The default value is 14 (4.17 rounded) per Palkovits, 1972.
%
%      GC_dendrite_length- This is the maximum length of a GC dendrite. 
%                  I have opted to use 20 um, because of the relative rarity of dendrites longer this.
%
%      GC_radius- This is the radius of an average GC.
%                  I have have opted to use 3 um. This is highly species dependent.
%
%      MFr_dist - This is the average distance between two neighboring MFrs, in um.
%                  This is 18.4um per Palkovits, 1972.
%
%%

GCL_Space = [100 100 250]; %We used 100 x 100 x 250 for most experiments. 
MFr_Density = 98800;
% Calculate MFr per space modeled.
MFr_Count = round(MFr_Density * (prod(GCL_Space)/1000^3)); %This needs to be a whole number.

GC_MFr_Convergence = 56;
GC_MFr_Divergence = 4;

GC_Count = round(MFr_Count * GC_MFr_Convergence/GC_MFr_Divergence); 
% We're going to use this one, because of how we distribute GCs below.

GC_dendrite_length = 20;
MFr_radius = 5;
GC_radius = 3;
% Calculate center to center distance, since we will be using point based modeling.
GC_MFr_centers_dist = GC_dendrite_length + MFr_radius + GC_radius; 
GC_MFr_c2c = MFr_radius + GC_radius;

MFr_dist = 18.4;

%% This set of functions generates the GCL system.
% GC and MFr are object-based lists of individual GCs and MFrs, with their coordinates,
% their 'identity' and their connectivity.
Relaxation = 2; % This is the allowed variance around MFr_dist, make it reasonable, 2 is default.
MFr = MFr_Generation(GCL_Space, MFr_Count, MFr_dist, Relaxation);

% Change this to 1 to see placement of MFrs.
Display = 0; 
if Display == 1;
    figure(); hold on; title(['C2C Dist is ',num2str(MFr_dist)]); daspect([1 1 1]);
    for imf = 1:size(MFr,2)
        scatter3(MFr(imf).x,MFr(imf).y,MFr(imf).z)
    end
end

Save = 0; % MFRs don't necessarily need to be saved, since they are included in the GrC system.
if Save == 1;
    descrip = sprintf('Saved_MFr_Systems/system of %.0i MFrs in %.0i by %.0i by %.0i on %s ex %i',...
                        MFr_Count, GCL_Space(1), GCL_Space(2), GCL_Space(3), date, randi([0 100000],1,1))
    mkdir('Saved_MFr_Systems')
    save(descrip,'MFr')
end


[GC,MFr,~] = GC_Generation(MFr, GC_Count, GC_MFr_centers_dist, GC_MFr_c2c, GC_MFr_Divergence, GC_MFr_Convergence);

% Change this to 1 to see placement of GCs with the MFrs.
Display = 0; 
if Display == 1;
    figure(); hold on; daspect([1 1 1]);
    for imf = 1:size(MFr,2)
        scatter3(MFr(imf).x,MFr(imf).y,MFr(imf).z,30,'green','filled')
    end
    for imf = 1:size(GC,2)
        scatter3(GC(imf).x,GC(imf).y,GC(imf).z,15,'red','filled')
    end
end

% These should be saved.
Save = 1;
if Save == 1;
    GC_MFr = {GC MFr}
    descrip = sprintf('Saved_GCL_Systems/system of %.0i MFrs with %.0i GCs in %.0i by %.0i by %.0i on %s ex %i',...
                        MFr_Count, GC_Count, GCL_Space(1), GCL_Space(2), GCL_Space(3), date, randi([0 100000],1,1))
    mkdir('Saved_GCL_Systems')
    save(descrip,'GC_MFr')
end

end
