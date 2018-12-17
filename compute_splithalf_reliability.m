function stimTraces = compute_splithalf_reliability(stimTraces)

rois = {'V1', 'V2', 'V3', 'V4'};

if ieNotDefined('stimTraces')
  saveName = '~/proj/texER/stimTraces.mat';
  computeStimTraces(rois, saveName);
  stimTraces = load(saveName);
end

%% Compute and plot correlation matrix.
conds = combvec(1:length(stimTraces.stimNames), 1:length(stimTraces.stimNames));

nObservations = [2,4,6,8,11];
corr_mtx = zeros(length(nObservations), length(rois), length(stimTraces.stimNames), length(stimTraces.stimNames));

for oi = 1:length(nObservations);
  nSmps = nObservations(oi);
  disp(sprintf('Computing correlation matrix with nObservations = %i', nSmps*2)); 
  rois = {'V1', 'V2', 'V3', 'V4'};
  for ri = 1:length(rois)
    traces = stimTraces.(sprintf('%s_traces', rois{ri}));

    disppercent(-inf, sprintf('ROI: %s', rois{ri}));
    for ci = 1:size(conds,1)
      cond = conds(ci,:);

      t1 = squeeze(mean(traces{cond(1)},2));
      t2 = squeeze(mean(traces{cond(2)},2));

      % Draw 1000 samples of prespecified size and take correlations.
      corrs = [];
      for i = 1:100
        idx = randperm(min(size(t2,2), size(t1,2))); % randomly permute the number of observations

        a1 = mean(t1(:, idx(1:nSmps)),2); % take nSmps observations of condition 1
        a2 = mean(t2(:, idx(nSmps+1 : 2*nSmps)), 2); % and a different nSmps observations of condition 2

        c = corrcoef(a1, a2); % and correlate them.
        corrs(i) = c(1,2);
      end
      corr_mtx(oi, ri, cond(1), cond(2)) = mean(corrs);
      
      disppercent(ci / size(conds,1));
    end
    disppercent(inf);

  end

  % Plot correlation matrix.
  plotCorMtx = 0;
  if plotCorMtx == 1
    h=figure; 
    for i = 1:4
      subplot(2,2,i);
      imagesc(stimTraces.(rois{i}).(sprintf('corr_mtx_%i', 2*nSmps)));
      set(gca, 'XTick', 1:21); set(gca, 'YTick', 1:21);
      set(gca, 'YTickLabel', stimTraces.stimNames); 
      colormap('hot');
      colorbar;
      title(sprintf('%s Correlation Matrix, nObs = %g', rois{i}, nSmps*2), 'FontSize', 18, 'FontName', 'Gill Sans');
    end
  end
  
end

%% Extract within and between family correlations
nConditions = length(stimTraces.stimNames);
cond_pairs = combvec(1:nConditions, 1:nConditions);
corr_cond = zeros(length(nObservations), length(rois), 4);
for oi = 1:length(nObservations)
  kObs = nObservations(oi);

  for ri = 1:length(rois)
    % Get corr_mtx
    cm = squeeze(corr_mtx(oi,ri,:,:));
    
    cond1 = []; cond2 = []; cond3 = []; cond4 = [];
    for ci = 1:size(cond_pairs,1)
      conds = cond_pairs(ci,:);

      fam1 = stimTraces.stimValues(1, conds(1)); fam2 = stimTraces.stimValues(1, conds(2));
      smp1 = stimTraces.stimValues(2, conds(1)); smp2 = stimTraces.stimValues(2, conds(2));

      if fam1 == fam2 % Within-Family
        if smp1 < 5 && smp2 < 5 % tex-tex
          cond1 = [cond1 squeeze(cm(conds(1), conds(2), :))];

        elseif xor(smp1<5, smp2<5) % tex-noise.
          cond2 = [cond2 squeeze(cm(conds(1), conds(2), :))];
        end

      else % Between Family
        if smp1 < 5 && smp2 < 5 % tex-tex
          cond3 = [cond3 squeeze(cm(conds(1), conds(2)))];
        elseif xor(smp1<5, smp2<5) % tex-noise.
          cond4 = [cond4 squeeze(cm(conds(1), conds(2)))];
        end
      end
    end

    for i = 1:4
      corr_cond(oi,ri,i) = mean(eval(sprintf('cond%i', i)));
    end
  end
end

%%  Get confidence intervals and chance reliability level.
confidence_interval = load('~/proj/texER/confidence_interval.mat');
mean_chance_corr = compute_splithalf_chance(stimTraces);

%% Plot within and between family split half correlations as a function of number of repeats.
 
colors = brewermap(5, 'YlGnBu');
colors = colors(2:end, :);

f = figure; set(gcf, 'Color', [1 1 1]);
%set(gcf, 'Position', [78.019, 4.015, 27.8967, 32.7928]);

for i = 1:4
  
  subplot(2,2,i);
  Y = corr_cond(:,:,i);
  ci = squeeze(confidence_interval.conf_intrvl(:,:,i,:));
  
  for j = 1:4
    myerrorbar(2*nObservations, Y(:,j), 'yLow', ci(:,j,1), 'yHigh', ci(:,j,2), 'Symbol=--o', 'Color', colors(j,:)); hold on;
  end
  ylim([-.1 .3]);
  title(sprintf('%s split-half correlation', confidence_interval.conditions{i}));
  xlabel('Number of repeated presentations');
  ylabel('Correlation');
  %if i==1; legend(rois);end
  
   hold on;
   for j = 1:4
     myerrorbar(2*nObservations, mean_chance_corr(:,j), 'Symbol=:.', 'Color', colors(j,:)); hold on;
   end
  drawPublishAxis;
end
%savepdf(f, '~/proj/texER/splithalf_corr.pdf')
%%
keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%
function cv = combvec(A,B)
[a,b] = meshgrid(A,B);
c = cat(2,a',b');
cv = reshape(c, [], 2);

%%%%%%%%%%%%%%%%%%%%%%%%%
% function to parse name
%%%%%%%%%%%%%%%%%%%%%%%%%
function [stimValues, condNames] = parseConditionName(stimNames)

if strcmp(class(stimNames), 'char')
  stimNames = {stimNames};
end

nConds = length(strsplit(stimNames{1}, ' and '));
stimValues = zeros(nConds, length(stimNames));
condNames = {};

for ci = 1:length(stimNames)
  conds = strsplit(stimNames{ci}, ' and ');
  for i = 1:nConds
    splits = strsplit(conds{i}, '=');
    stimValues(i,ci) = str2num(splits{2});
    condNames{i} = splits{1};
  end
end
