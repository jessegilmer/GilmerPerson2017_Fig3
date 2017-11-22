%% Created 11/23/16 by Gilmer, J.
% This is the MFr generator. It needs to be fed:
%                          Space - The modelling space that should be defined in GCL_Main.
%
%                          Count - The number of MFrs to create.
%
%                          MFr_dist- The mean distance between MFrs, as defined in GCL_Main.
%
%                          Relaxation - How much you want the MFr_dist requirement to be relaxed by (in um)
%
%   Example: MF = MFr_Generation([1000 1000 1000], 90000, 18.4, 2)

%%
function [MFr] = MFr_Generation(Space,Count,Center_dist,Relaxation)
tic;
rng('shuffle');
metatime = 0

% These are the empty arguement catches. Please just use 4 arguements.
if (nargin == 0)
    Space = [100,100,100];
    % You'll probably want to change this to what you want.
    Count = round(98800 * (prod(Space)/1000^3));
    Relaxation = .5;
    Center_dist = 18.4;
    
elseif (nargin == 1)
    Count = round(98800 * (prod(Space)/1000^3));
    Relaxation = .5;
    Center_dist = 18.4;
    
elseif (nargin == 2)
    Center_dist = 18.4;
    Relaxation = .5;
    
elseif (nargin == 3)
    Relaxation = .5;
    
end


%% Initial Point Generation

% Set boundary checks.
xlim = Space(1);
ylim = Space(2);
zlim = Space(3);

%Protect against boundary issues
div = 10;
div_errorx = xlim - (xlim - xlim/div);
div_errory = ylim - (ylim - ylim/div);
div_errorz = zlim - (zlim - zlim/div);
div_error = min([div_errorx div_errory div_errorz]);
while div_error < Center_dist
    div = div - 1;
    div_errorx = xlim - (xlim - xlim/div);
    div_errory = ylim - (ylim - ylim/div);
    div_errorz = zlim - (zlim - zlim/div);
    div_error = min([div_errorx div_errory div_errorz]);
end

% MFr(1) is going to be placed in the center of the given space, with some variance.
MFr(1).x = ceil(Space(1)/2) + randi([ceil(-1000*(xlim/div)) floor(1000*(xlim/div))],1,1)/1000;
MFr(1).y = ceil(Space(2)/2) + randi([ceil(-1000*(ylim/div)) floor(1000*(ylim/div))],1,1)/1000;
MFr(1).z = ceil(Space(3)/2) + randi([ceil(-1000*(zlim/div)) floor(1000*(zlim/div))],1,1)/1000;



MFr(1).id = 1;

%daughters will be used when GCs are generated.
MFr(1).daughters = 0;
MFr(1).daughterid = [];
MFr(1).neighbors = 1;

% MFr(2) is going to be Center_dist away from MF(1) at some random set of
% angles, with some variance.
p = Center_dist +  randi([ceil(-1000*Relaxation) floor(1000*(Relaxation))],1,1)/1000;
phi = 2*pi*(randi([0 1000],1,1)/1000);
ang = 2*pi*(randi([0 1000],1,1)/1000);

MFr(2).x = MFr(1).x + p*sin(phi)*cos(ang);
MFr(2).y = MFr(1).y + p*sin(phi)*sin(ang);
MFr(2).z = MFr(1).z + p*cos(phi);

MFr(2).id = 2;
MFr(2).daughters = 0;
MFr(2).daughterid = [];
MFr(2).neighbors = 1;
clear a  b c phi ang;

%set a variable to keep track of MFrs created.
mfct = 2;
mfr_n_list = zeros(1,Count);
mfr_n_list(1) = 1;
mfr_n_list(2) = 1;
runners(1,:) = [MFr(1).x MFr(1).y MFr(1).z];
runners(2,:) = [MFr(2).x MFr(2).y MFr(2).z];

%% Generative process for MFr(3) to MFr(Count)
while mfct < Count;
    
    % Choose a random MFr to use as a point reference, and place it
    % Center_dist +- variance away at some random set of angles.
    
    neig_check = mfr_n_list(1:mfct);
    neig_check = find(neig_check < 12); %theorhetically this is up to 18, but use like... 12 or so;
    randmf = neig_check(randi([1 length(neig_check)],1,1));
    
    p = Center_dist + randi([-Relaxation*1000 Relaxation*1000],1,1)/1000;
    phi = 2*pi*(randi([0 1000],1,1)/1000);
    ang = 2*pi*(randi([0 1000],1,1)/1000);
    x = MFr(randmf).x + p*sin(phi)*cos(ang);
    y = MFr(randmf).y + p*sin(phi)*sin(ang);
    z = MFr(randmf).z + p*cos(phi);
    
    metadist = (sqrt((x - runners(:,1)).^2 + (y - runners(:,2)).^2 + (z - runners(:,3)).^2));
    clear alpha
    alpha = find(metadist < Center_dist + Relaxation);
    metadist = sort(metadist);
    if mfct > 2 %This needs at least 2 MFrs
        while min(metadist) < Center_dist - Relaxation || metadist(2) > Center_dist + Relaxation || metadist(3) > Center_dist + Relaxation
            randmf = neig_check(randi([1 length(neig_check)],1,1));
            p = Center_dist + randi([-Relaxation*1000 Relaxation*1000],1,1)/1000;
            phi = 2*pi*(randi([0 1000],1,1)/1000);
            ang = 2*pi*(randi([0 1000],1,1)/1000);
            x = MFr(randmf).x + p*sin(phi)*cos(ang);
            y = MFr(randmf).y + p*sin(phi)*sin(ang);
            z = MFr(randmf).z + p*cos(phi);
            
            metadist = (sqrt((x - runners(:,1)).^2 + (y - runners(:,2)).^2 + (z - runners(:,3)).^2));
            clear alpha
            alpha = find(metadist < Center_dist + Relaxation);
            metadist = sort(metadist);
        end
    else
        while min(metadist) < Center_dist - Relaxation || metadist(2) > Center_dist + Relaxation
            randmf = neig_check(randi([1 length(neig_check)],1,1));
            p = Center_dist + randi([-Relaxation*1000 Relaxation*1000],1,1)/1000;
            phi = 2*pi*(randi([0 1000],1,1)/1000);
            ang = 2*pi*(randi([0 1000],1,1)/1000);
            x = MFr(randmf).x + p*sin(phi)*cos(ang);
            y = MFr(randmf).y + p*sin(phi)*sin(ang);
            z = MFr(randmf).z + p*cos(phi);
            
            metadist = (sqrt((x - runners(:,1)).^2 + (y - runners(:,2)).^2 + (z - runners(:,3)).^2));
            clear alpha
            alpha = find(metadist < Center_dist + Relaxation);
            metadist = sort(metadist);
        end
    end
    
    
    % Make sure this point is within Space;
    if x < 0.1 || x > xlim
        continue;
    end
    if y < 0.1 || y > ylim
        continue;
    end
    if z < 0.1 || z > zlim
        continue;
    end
    
    
    mfct = mfct + 1;
    
    % To keep track of progress.
    if mod(mfct,1000) == 0;
        progress = [mfct Count]
        time = toc
        metatime = metatime + time
        tic
    end
    
    %store good results in the MFr object.
    MFr(mfct).x = x;
    MFr(mfct).y = y;
    MFr(mfct).z = z;
    MFr(mfct).id = mfct;
    MFr(mfct).daughters = 0;
    MFr(mfct).neighbors = length(alpha);
    runners(mfct,:) = [x y z];
    
end
end
