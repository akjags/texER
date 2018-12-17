function analyzeTextureTraces()

%% Set up view
v = newView;
v = viewSet(v, 'curGroup', 'Concatenation');
v = viewSet(v, 'curScan', 1);

%% Load analysis and get d variable
v = loadAnalysis(v, 'erAnal/erAnal2.mat');

d = viewGet(v, 'd');
overlays = viewGet(v, 'overlays');
r2 = overlays(1).data{1}; % Get r2 map

% Get traces and condition names
ehdr = d.ehdr;
stimNames = d.stimNames;
[stimValues, condNames] = parseConditionName(stimNames);

cutoff=.1;
%% Loop through ROIs
rois = {'V1', 'V2', 'V3', 'V4'};
keyboard
max_tex = []; max_noise = [];
for ri = 1:length(rois)
  roi = loadROITSeries(v, rois{ri});
  
  tex_trc = []; noise_trc = [];
  trcs = []; modIdx = [];
  for vi = 1:size(roi.scanCoords,2)
    sc = roi.scanCoords(:,vi);
    
    if r2(sc(1), sc(2), sc(3)) > cutoff
      trc = squeeze(ehdr(sc(1), sc(2), sc(3), :, :));
      
      texAmp = max(trc([2 4 6],:), [], 2);
      noiseAmp = max(trc([1 3 5],:), [], 2);
      modIdx = [modIdx, mean(texAmp - noiseAmp)./mean(texAmp+noiseAmp)];
      trcs= cat(3,trc, trcs);
    end
  end
  tex_trc = trcs([2 4 6],:,:);
  noise_trc = trcs([1 3 5],:,:);
  
  mean_tex = mean(tex_trc, 3);
  mean_noise = mean(noise_trc,3);
  
  max_tex(ri) = mean(mean(max(mean_tex, [], 2)));
  max_noise(ri) = mean(mean(max(mean_noise, [], 2)));
  
  disp(sprintf('%s Texture Amplitude: %g, Noise Amplitude %g', roi.name, max_tex(ri), max_noise(ri)));
  mod_idxs(ri) = mean(modIdx);
end
  
%% Plot amplitudes
amps = [max_tex; max_noise];
mod_idx2 = (amps(1,:) - amps(2,:)) ./ (amps(1,:) + amps(2,:));


figure; set(gcf, 'Color', [1 1 1]);
subplot(3,1,1);
bar(amps'); colormap('Winter');
title('Amplitude of Average Traces', 'FontName', 'Gill Sans', 'FontSize', 18);
set(gca, 'XTickLabel', rois);
legend({'Texture', 'Noise'});
box off;

subplot(3,1,2);
bar(mod_idxs); colormap('Winter'); box off;
title('Modulation Index (averaged across voxels)', 'FontName', 'Gill Sans', 'FontSize', 18);
set(gca, 'XTickLabel', rois);

subplot(3,1,3);
bar(mod_idx2); colormap('Winter'); box off;
title('Modulation Index of Average Trace', 'FontName', 'Gill Sans', 'FontSize', 18);
set(gca, 'XTickLabel', rois);

%%
% returns the k largest value in an array A, along dimension dim.
function kLargest = maxk(A, k, dim)

sorted = sort(A,dim, 'descend');

kLargest = sorted(:, 1:k);

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
