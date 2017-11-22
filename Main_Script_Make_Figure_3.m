%% Guide to making Figure 3 of Gilmer and Person, 2017.
% Written by J Gilmer. 11/21/17
% Enjoy.

% This clears the workspace.
clear all; close all;

% First, you will need a set of simulations. 
% This will probably take quite a while.

n_systems = 150; % Set the number of systems you want to generate. 
                 % 150 is typically what we used.
                 
for i = 1:n_systems
    
    GCL_Main();
    % GCL_Main has the parameters of the physiology.
    % It calls MFr_Generation to populate the MFRs.
    % It then calls GC_Generation to populate the GrCs.
    % (Sorry about the shifting terminology.)
    % You'll also get a readout of the progress in the command window.
    
end
% You should now have a folder full of generated systems.

%% Now we want to generate the CE metric

% Reclear the workspace:
clear all; close all;

step_size = 1; % This is the change in MFR diversity for the function.
trials = 5; % This is how many repeats of the function to run per system.
% (if step_size isn't 1 you're going to have to do some data allignment when
% we process the dimensionality. i.e you need to specify what MFR diversity
% levels were generated in the data processing loop.)

GCdirect = dir('Saved_GCL_Systems'); % Target the saved systems folder.

for i = 3:length(GCdirect) % The first two paths are always just some dots.
    
    load(sprintf('%s/%s',GCdirect(i).folder,GCdirect(i).name)) % Load in a system.
    GC = GC_MFr{1}; %Reassign the GrCs and MFRs;
    MFr = GC_MFr{2};
    clear GC_MFr;
    
    for j = 1:5
        % You should probably set this to repeat a few times, so you resample the same system.
        
        % It's also probably worth noting that this processes takes a 
        % LONG time to run...
        
        dim_result = Codon_Dim(MFr,GC,trials,step_size,0);
        % This prints a progress log, of [trialnumber and MFRdiversity] for both 'conditions'.
        
        descrip = sprintf('Dim Data/Codon_Dim_Result session %i trial %i',i,j)
        mkdir('Dim Data')
        save(descrip)
    end
    
    close all;
    clearvars -except GCdirect i trials step_size;
end

%% Now we need to extract the data:

% Reclear the workspace:
clear all; close all;

dimdirect = dir('Dim Data'); % Target the dimensionality folder.


for ii = 3:length(dimdirect) % The first two paths are always just some dots.
    
    load(sprintf('%s/%s',dimdirect(ii).folder,dimdirect(ii).name))% Load in data
    
    for jj = 1:length(dim_result)
        
        trial = dim_result(jj).trial;
        
        if ~isempty(trial) %This should give you a little flexibility...
            
            for k = 1:length(trial)
                dspatial(k) = trial(k).dim; %this is the spatial dimensionality.
                dnonspatial(k) = trial(k).marrdim; % this is the non-spatial dimensionality.
            end
            
            % We want the mean dimensionality per system:
            dimension_spatial(ii-2,jj) = mean(dspatial); 
            dimension_nonspatial(ii-2,jj) = mean(dnonspatial);
            clear dspatial dnonspatial;
        end
        
    end
    
    clearvars -except dimdirect i dimension_spatial dimension_nonspatial;
end

% Finally, we get the CE metric:
dim_spatial = mean(dimension_spatial);
dim_nonspatial = mean(dimension_nonspatial);

figure(1);
hold on;
plot(dim_spatial,'r');
plot(dim_nonspatial,'k');
legend('spatial dimensionality','non-spatial dimensionality')
xlabel('MFR diversity')
ylabel('CE Dimensionality')

%% Now for TECE:

% Reclear the workspace:
clearvars -except dim_spatial dim_nonspatial; 

time = 500; % PKJ integration time.
unit = 20; % GrC response time.
q = 1:100; % MFR diversity level.
M = 3458; %you'll want to change this if you have more or less GrCs.
bestR = (time/unit); %This is the best potential redundancy.

dim = dim_spatial;
dimM = dim_nonspatial;

split = M./q; %division of GrCs per MFR.
maxbest = M ./ bestR;

dimr = M./dim;

% To determine waste:
% dim.*bestR is the number of unique units times the nubmer of
% explansions that are optimal. Those will produce diversity,
% everything else is waste. therefore R - that number.
r_waste = M - (time/unit .* dim);
r_waste(r_waste <= 0) = 0;

r_loss = ((time/unit) - M./dim) .* M*unit/time;
r_loss(r_loss < 0) = 0;


lw = M - (r_loss + r_waste);

maxbest = M ./ bestR;

dimr = M./dimM;

% To determine waste:
% dim.*bestR is the number of unique units times the nubmer of
% explansions that are optimal. Those will produce diversity,
% everything else is waste. therefore R - that number.
r_waste = M - (time/unit .* dimM);
r_waste(r_waste <= 0) = 0;

r_loss = ((time/unit) - M./dimM) .* M*unit/time;
r_loss(r_loss < 0) = 0;

lw2 = M - (r_loss + r_waste);

plotting = 0; % You probably don't care what the losswaste function looks like.
if plotting == 1;
figure(3);
hold on;
plot(lw);
plot(lw2);
end

TECE1 = lw;
TECE2 = lw2;
figure(2);
hold on;
plot(TECE1, 'r')
plot(TECE2, 'k')
legend('spatial TECE dimensionality','non-spatial TECE dimensionality')
xlabel('MFR diversity')
ylabel('TECE Total Dimensionality')






