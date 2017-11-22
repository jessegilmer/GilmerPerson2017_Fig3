%Readme
%{
Created 6/29/17 by Gilmer, J.
This examines variance and correlation in granule cells.
Example: [Results] = Codon_Correlation(MFr,GC,1,1)
%}

function [ubervect] = Codon_Dim(MFr,GC,trials,step_size,Display)
% This thing will store trial information.
ubervect = struct();

% Save base states for re-loading in script.
MFr2 = MFr;
GC2 = GC;

% Plot display toggle.
if Display == 1;
    figure();
    hold on;
end

%% Limit generation size to 100.
if size(MFr,2) > 100;
    metasize = 100;
else
    metasize = size(MFr,2);
end

%% Mainscript loop.
for thr = 3
    for jj = 1:trials
        %% Control conditions loop.
        for ii = 1:step_size:100;
            [jj ii]
            idlow = 1;
            idhigh = ii;
            
            % make a randomized pool of MFr identities.
            pool_len = 0;
            offset = 0;
            while pool_len ~= ii;
                identity_pool = randi([idlow idhigh+offset],1,length(MFr2));
                pool_len = length(unique(identity_pool));
                if pool_len > ii
                    offset = offset - 1;
                else
                    offset = offset + 1;
                end
            end
            
            % Reload initial system.
            MFr = MFr2;
            GC = GC2;
            
            % Resassign codons.
            for i = 1:length(MFr);
                MFr(i).id2 = identity_pool(i);
                daughters = MFr(i).daughterid;
                for j = 1:length(daughters);
                    parents = GC(daughters(j)).parents;
                    for kk = 1:length(parents);
                        if parents(kk) == MFr(i).id;
                            parents(kk) = MFr(i).id2;
                        end
                    end
                    GC(daughters(j)).parents = parents;
                end
            end
            
            for i = 1:length(GC);
                codon(i,:) = GC(i).parents;
            end
            
            tic;
            % give random activity for 1000 points;
            act_log = zeros(1000,length(GC));
            activekey = 1:length(unique(identity_pool));
            xx = .9;            
            for t = 1:1000;
                actives = rand(1,length(unique(identity_pool)));
                actives = actives >= xx;
                gc_act = codon;
                for j = activekey
                    gc_act(gc_act == j) = actives(j);
                end
                act_log(t,:) = sum(gc_act');
                clear gc_act;
            end
            %threshhold is 3
            act_log(act_log<thr) = 0;
            act_log(act_log>=thr) = 1;
            
            corr_log = corrcoef(act_log);
            
            
            for i = 1:length(GC);
                codon_in(i,:) = GC(i).parents;
            end
            
            cod1_hold = codon_in;
            codon_in = unique(codon_in,'rows');
            %
            ct3 = 1;
            
            for i = 1:size(codon_in,1);
                cod_in = codon_in(i,:);
                sub_selector = nchoosek(cod_in,3);
                for j = 1:length(sub_selector);
                    codon3_in(ct3,:) = sub_selector(j,:);
                    ct3 = ct3 + 1;
                end
            end
            
            ct2 = 1;
            
            for i = 1:size(codon_in,1);
                cod_in = codon_in(i,:);
                sub_selector2 = nchoosek(cod_in,2);
                for j = 1:length(sub_selector2);
                    codon2_in(ct2,:) = sub_selector2(j,:);
                    ct2 = ct2 + 1;
                end
            end
            
            ct1 = 1;
            
            for i = 1:size(codon_in,1);
                cod_in = codon_in(i,:);
                sub_selector3 = nchoosek(cod_in,1);
                for j = 1:length(sub_selector3);
                    codon1_in(ct1,:) = sub_selector3(j,:);
                    ct1 = ct1 + 1;
                end
            end
            
            codon_in = cod1_hold;
            codon_in = sort(codon_in,2);
            c4_unique = unique(codon_in,'rows');
            
            codon2_in = sort(codon2_in,2);
            c2_unique = unique(codon2_in,'rows');
            
            codon3_in = sort(codon3_in,2);
            c3_unique = unique(codon3_in,'rows');
            
            c1_unique = unique(codon1_in);
            
            for jjj = 1:size(c4_unique,1)
                codon4_ck = c4_unique(jjj,:);
                out4_check(jjj) = sum(sum((codon_in(:,:) == codon4_ck(1,:))') == 4);
                cod4_registry(jjj,:) = codon4_ck;
            end
            
            for jjj = 1:size(c3_unique,1)
                codon3_ck = c3_unique(jjj,:);
                out3_check(jjj) = sum(sum((codon_in(:,:) == [codon3_ck(1,:) 0])') >= 3);
                cod3_registry(jjj,:) = codon3_ck;
            end
            
            for jjj = 1:size(c2_unique,1)
                codon2_ck = c2_unique(jjj,:);
                out2_check(jjj) = sum(sum((codon_in(:,:) == [codon2_ck(1,:) 0 0])') >= 2);
                cod2_registry(jjj,:) = codon2_ck;
            end
            
            for jjj = 1:size(c1_unique,1)
                codon1_ck = c1_unique(jjj,:);
                out1_check(jjj) = sum(sum((codon_in(:,:) == [codon1_ck(1,:) 0 0 0]) >= 1));
                cod1_registry(jjj,:) = codon1_ck;
            end
            
            ubervect(ii).trial(jj).out1 = out1_check;
            ubervect(ii).trial(jj).out2 = out2_check;
            ubervect(ii).trial(jj).out3 = out3_check;
            ubervect(ii).trial(jj).out4 = out4_check;
            
            ubervect(ii).trial(jj).reg1 = cod1_registry;
            ubervect(ii).trial(jj).reg2 = cod2_registry;
            ubervect(ii).trial(jj).reg3 = cod3_registry;
            ubervect(ii).trial(jj).reg4 = cod4_registry;
            
            toc
            
            corr_log(diag(1:length(act_log))>0) = NaN;
            
            ubervect(ii).trial(jj).actmean = nanmean(nanmean(act_log));
            ubervect(ii).trial(jj).corr = nanmean(nanmean(corr_log));
            ubervect(ii).trial(jj).var = nanmean(nanvar(corr_log));
            %         ubervect(ii).trial(jj).var = 1/ii;
            ubervect(ii).trial(jj).dim = ...
                1/((1/length(GC))+ubervect(ii).trial(jj).corr^2 + ubervect(ii).trial(jj).var);
            clearvars -except thr metasize Display step_size trials ii jj MFr2 GC2 ubervect;
            
        end
        %% Marr conditions loop.
        for ii = 1:step_size:100;
            [jj ii]
            idlow = 1;
            idhigh = ii;
            
            % Reload initial system.
            MFr = MFr2;
            GC = GC2;
            
            % Resassign codons.
            for i = 1:length(GC);
                GC(i).parents = randi([1 ii],1,4);
            end
            
            for i = 1:length(GC);
                codon(i,:) = GC(i).parents;
            end
            tic;
            
            % give random activity for 100 points;
            act_log = zeros(1000,length(GC));
            activekey = 1:ii;
            xx = .9;
            for t = 1:1000;
                actives = rand(1,ii);
                actives = actives >= xx;
                gc_act = codon;
                for j = activekey
                    gc_act(gc_act == j) = actives(j);
                end
                act_log(t,:) = sum(gc_act');
                clear gc_act;
            end
            %threshhold is 3
            act_log(act_log<thr) = 0;
            act_log(act_log>=thr) = 1;
            
            corr_log = corrcoef(act_log);
            
            for i = 1:length(GC);
                codon_in(i,:) = GC(i).parents;
            end
            
            cod1_hold = codon_in;
            codon_in = unique(codon_in,'rows');
            %
            ct3 = 1;
            
            for i = 1:size(codon_in,1);
                cod_in = codon_in(i,:);
                sub_selector = nchoosek(cod_in,3);
                for j = 1:length(sub_selector);
                    codon3_in(ct3,:) = sub_selector(j,:);
                    ct3 = ct3 + 1;
                end
            end
            
            ct2 = 1;
            
            for i = 1:size(codon_in,1);
                cod_in = codon_in(i,:);
                sub_selector2 = nchoosek(cod_in,2);
                for j = 1:length(sub_selector2);
                    codon2_in(ct2,:) = sub_selector2(j,:);
                    ct2 = ct2 + 1;
                end
            end
            
            ct1 = 1;
            
            for i = 1:size(codon_in,1);
                cod_in = codon_in(i,:);
                sub_selector3 = nchoosek(cod_in,1);
                for j = 1:length(sub_selector3);
                    codon1_in(ct1,:) = sub_selector3(j,:);
                    ct1 = ct1 + 1;
                end
            end
            
            codon_in = cod1_hold;
            codon_in = sort(codon_in,2);
            c4_unique = unique(codon_in,'rows');
            
            codon2_in = sort(codon2_in,2);
            c2_unique = unique(codon2_in,'rows');
            
            codon3_in = sort(codon3_in,2);
            c3_unique = unique(codon3_in,'rows');
            
            c1_unique = unique(codon1_in);
            
            for jjj = 1:size(c4_unique,1)
                codon4_ck = c4_unique(jjj,:);
                out4_check(jjj) = sum(sum((codon_in(:,:) == codon4_ck(1,:))') == 4);
                cod4_registry(jjj,:) = codon4_ck;
            end
            
            for jjj = 1:size(c3_unique,1)
                codon3_ck = c3_unique(jjj,:);
                out3_check(jjj) = sum(sum((codon_in(:,:) == [codon3_ck(1,:) 0])') >= 3);
                cod3_registry(jjj,:) = codon3_ck;
            end
            
            for jjj = 1:size(c2_unique,1)
                codon2_ck = c2_unique(jjj,:);
                out2_check(jjj) = sum(sum((codon_in(:,:) == [codon2_ck(1,:) 0 0])') >= 2);
                cod2_registry(jjj,:) = codon2_ck;
            end
            
            for jjj = 1:size(c1_unique,1)
                codon1_ck = c1_unique(jjj,:);
                out1_check(jjj) = sum(sum((codon_in(:,:) == [codon1_ck(1,:) 0 0 0]) >= 1));
                cod1_registry(jjj,:) = codon1_ck;
            end
            
            ubervect(ii).trial(jj).out1marr = out1_check;
            ubervect(ii).trial(jj).out2marr = out2_check;
            ubervect(ii).trial(jj).out3marr = out3_check;
            ubervect(ii).trial(jj).out4marr = out4_check;
            
            ubervect(ii).trial(jj).reg1 = cod1_registry;
            ubervect(ii).trial(jj).reg2 = cod2_registry;
            ubervect(ii).trial(jj).reg3 = cod3_registry;
            ubervect(ii).trial(jj).reg4 = cod4_registry;
            
            
            toc
            
            corr_log(diag(1:length(act_log))>0) = NaN;
            
            ubervect(ii).trial(jj).marractmean = nanmean(nanmean(act_log));
            ubervect(ii).trial(jj).marrcorr = nanmean(nanmean(corr_log));
            ubervect(ii).trial(jj).marrvar = nanmean(nanvar(corr_log));
            ubervect(ii).trial(jj).marrdim = ...
                1/((1/length(GC))+ubervect(ii).trial(jj).marrcorr^2 + ubervect(ii).trial(jj).marrvar);
            clearvars -except thr metasize Display step_size trials ii jj MFr2 GC2 ubervect;
        end
    end
end
