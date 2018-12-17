function drawEventRelatedFigs(rois, saveFigs, r2cutoff)

if ieNotDefined('rois')
  rois = {'V1', 'V2', 'V3', 'V4'};
end
if ieNotDefined('saveFigs')
  saveFigs = 0;
end
if ieNotDefined('saveLoc')
  saveLoc = '~/proj/texER/Figures2'; 
end
if ieNotDefined('r2cutoff')
  r2cutoff = 0.2;
end

%% Load pre-computed stimulus traces.
saveName = '~/proj/texER/stimTraceAmps.mat';
if ~exist(saveName)
  computeStimTraces(rois, saveName);
end
stimTraces = load(saveName);
stimTraces2 = load('~/proj/texER/stimTraces.mat');
stimValues = stimTraces.stimValues; % 2x21 array
condNames = stimTraces.condNames; % Conditions: Texture Family x Sample Index
ehdr = stimTraces.ehdr; % hemodynamic response per condition
scanDims = size(ehdr);
scanDims = scanDims(1:3); 

texIdx = [1 2 3 4]; % which samples are textures (textures = 1-4, noise = 5-7)
isTexture = arrayfun(@(x) any(x==texIdx), stimValues(2,:));

texFams = unique(stimValues(1,:));


colors = brewermap(5, 'YlGnBu');
colors = colors(2:end, :);

keyboard

%% First plot average timecourse (ehdr)
iehdr = reshape(ehdr, prod(scanDims), size(ehdr, 4), size(ehdr, 5));

for ri = 1:length(rois)
  sc = stimTraces.(rois{ri}).scanCoords;
  
  % select only the scan coords whose r2 exceeds cutoff
  sc = sc(:, stimTraces.(rois{ri}).r2 > r2cutoff);
  ind = sub2ind(scanDims, sc(1,:), sc(2,:), sc(3,:)); % convert scan coords to linear index

  % Get all ehdr's and take mean across all voxels.
  roi_ehdr = iehdr(ind,:,:); % nVoxels x 21 x 50
  
  for tfI = 1:length(texFams)
    texTraces = squeeze(mean(roi_ehdr(:,stimValues(1,:) == texFams(tfI) & isTexture,:), 2));
    noiseTraces = squeeze(mean(roi_ehdr(:,stimValues(1,:) == texFams(tfI) & ~isTexture,:), 2));
    
    txTrc = mean(texTraces, 1);
    nsTrc = mean(noiseTraces, 1);
    count = (ri-1)*(length(texFams)+1)+tfI;
    
    subplot(length(rois), length(texFams)+1, count);
    h2=myerrorbar((1:length(nsTrc))/2, nsTrc, 'Color=b', 'yError', 1.96*std(noiseTraces,[],1)/sqrt(size(noiseTraces,1)), 'yErrorBarType=fill','fillAlpha=0.2');
    h1=myerrorbar((1:length(txTrc))/2, txTrc, 'Color=g', 'yError', 1.96*std(texTraces,[],1)/sqrt(size(texTraces,1)), 'yErrorBarType=fill','fillAlpha=0.2'); hold on;
    set(gca, 'XTick',1:5:(length(txTrc)/2));
    set(gca, 'XTickLabel', 1:5:(length(txTrc)/2));
    xlabel('Time (sec)'); ylabel('% signal change');
    ylim([-1 2]);
    title(sprintf('%s, %s', rois{ri}, stimTraces.imNames{tfI}));
    if count == 1
      legend([h1,h2],{'Textures', 'Phase-Scrambled'}, 'Location', 'Northwest');
    end
    drawPublishAxis('labelFontSize', 14);
  end
  % Plot mean across all texture families
  texTraces = squeeze(mean(roi_ehdr(:,isTexture,:), 2));
  noiseTraces = squeeze(mean(roi_ehdr(:,~isTexture,:), 2));

  txTrc = mean(texTraces, 1);
  nsTrc = mean(noiseTraces, 1);
  count = (ri-1)*(length(texFams)+1)+tfI;

  count = (ri-1)*(length(texFams)+1)+length(texFams)+1;
  subplot(length(rois), length(texFams)+1, count);
  myerrorbar((1:length(nsTrc))/2, nsTrc, 'Color=b', 'yError', 1.96*std(mean(noiseTraces, 1),[],3)/sqrt(size(noiseTraces,3)), 'yErrorBarType=fill','fillAlpha=0.2');
  myerrorbar((1:length(txTrc))/2, txTrc, 'Color=g', 'yError', 1.96*std(mean(texTraces, 1),[],3)/sqrt(size(texTraces,3)), 'yErrorBarType=fill','fillAlpha=0.2'); hold on;
  set(gca, 'XTick',1:5:(length(txTrc)/2));
  set(gca, 'XTickLabel', 1:5:(length(txTrc)/2));

  xlabel('Time (sec)');
  ylabel('% signal change');
  ylim([-1 2]);
  title(sprintf('%s - All', rois{ri}));
  drawPublishAxis('labelFontSize', 14);
end
set(gcf, 'Position', [67.6, 12.6, 45.8, 26.2]);

savepdf(f, [saveLoc '/timecourses']);
saveas(f, [saveLoc '/timecourses.png']);

% Plot mean across images
%%
keyboard
%% Second, plot the modulation index for each roi

[mod_index, mod_index_CI] = ER_computeModulationIndex(rois, stimTraces, r2cutoff);

f = figure;
% First plot as a bar graph averaged across all texture families
subplot(2,1,1);
avg_mod_idx = mean(mod_index, 2);
for i = 1:length(rois)
  bar(i, avg_mod_idx(i), 'FaceColor', colors(i,:)); hold on;
end
myerrorbar(1:length(rois), mean(mod_index,2), 'yError', mean(mod_index_CI, 2), 'Symbol', '.','MarkerEdgeColor', 'k', 'Color', 'k');
set(gca, 'XTick', 1:length(rois));
set(gca, 'XTickLabel', rois);
title('Modulation Index');
drawPublishAxis('labelFontSize', 14);

% Then split it out by texture family
subplot(2,1,2);
handles=[];
for i = 1:length(rois)
  h = myerrorbar(1:length(texFams), mod_index(i,:), 'yError', mod_index_CI(i,:), 'Symbol', 'o', 'Color', colors(i,:), 'MarkerSize', 10);
  handles = [handles h];
  hold on;
end
set(gca, 'XTick', 1:length(texFams));
set(gca, 'XTickLabel', stimTraces.imNames);
title('Modulation Index, split by texture family');
legend(handles, rois, 'Location', 'Northwest');
legend('boxoff');
drawPublishAxis('labelFontSize', 14);
set(f, 'Position', [0.11 16.23 21.2 21.66]);

if saveFigs
  savepdf(f, [saveLoc '/modulation_index.pdf']);
end
%%

%% Third, plot the within family texture correlation, by texture class.
corr_mtx = ER_getCorrMtx(rois, stimTraces, r2cutoff);
corrs = ER_getCorrStruct(rois, stimTraces, corr_mtx);

f=figure;
% First plot bar plots averaged across texture families
subplot(2,1,1);
col = colors(3,:);
bar((1:length(rois))-.125, mean(corrs.within_textex,2), .25, 'FaceColor', col);  hold on;
bar((1:length(rois))+.125, mean(corrs.between_textex,2), .25, 'FaceColor', col+(1-col)*(.5)); 
myerrorbar((1:length(rois))-.125, mean(corrs.within_textex,2), 'yError', mean(corrs.within_textex_se,2), 'Symbol', '.');
myerrorbar((1:length(rois))+.125, mean(corrs.between_textex,2), 'yError', mean(corrs.between_textex_se,2), 'Symbol', '.');
ylim([-.1, .3]);
set(gca, 'XTick', 1:length(rois));
set(gca, 'XTickLabel', rois);
title('Within vs Between family texture correlation');
legend({'Within', 'Between'}, 'Location', 'Northwest');
legend('boxoff');
drawPublishAxis('labelFontSize', 14);

% Next plot point plots split out by texture family.
subplot(2,1,2);
handles = [];
for ti = 1:length(texFams)
  col = colors(ti,:);
  h = myerrorbar((1:length(rois))-0.1, corrs.('within_textex')(:,ti), 'yError', corrs.('within_textex_se')(:,ti), 'Symbol', 'o', 'Color', col, 'MarkerSize', 10); hold on;
  handles = [handles h];
  h = myerrorbar((1:length(rois))+0.1, corrs.('between_textex')(:,ti), 'yError', corrs.('between_textex_se')(:,ti), 'Symbol', 'o', 'Color', col+(1-col)*(.5), 'MarkerSize', 10); hold on;
  handles = [handles h];
end
ylim([-.1, .3]);
set(gca, 'XTick', 1:length(rois));
set(gca, 'XTickLabel', rois);
title('Within vs Between family texture correlation, split by texture family');
xlim([0 length(rois)+1]);
labels=[cellfun(@(x) [x ' Within'], stimTraces.imNames, 'un', 0), cellfun(@(x) [x ' Between'], stimTraces.imNames, 'un', 0)];
legend(handles([1:2:end, 2:2:end]), labels, 'Location', 'Northwest');
legend('boxoff');
drawPublishAxis('labelFontSize', 14);
set(f, 'Position', [21.4, 16.2, 23.1, 21.2]);
if saveFigs
  savepdf(f, [saveLoc '/withinbetween_texture_correlation.pdf']);
end
%% Fourth, plot within family texture-noise correlation, by texture class.
f=figure;
% First plot bar plots averaged across texture families
subplot(2,1,1);
col = colors(3,:);
bar((1:length(rois))-.125, mean(corrs.within_texnoise,2), .25, 'FaceColor', col);  hold on;
bar((1:length(rois))+.125, mean(corrs.between_texnoise,2), .25, 'FaceColor', col+(1-col)*(.5)); 
myerrorbar((1:length(rois))-.125, mean(corrs.within_texnoise,2), 'yError', mean(corrs.within_textex_se,2), 'Symbol', '.');
myerrorbar((1:length(rois))+.125, mean(corrs.between_texnoise,2), 'yError', mean(corrs.between_textex_se,2), 'Symbol', '.');
ylim([-.1, .3]);
set(gca, 'XTick', 1:length(rois));
set(gca, 'XTickLabel', rois);
title('Correlation between textures and phase-scrambles (within vs between family)');
legend({'Within', 'Between'}, 'Location', 'Northwest');
legend('boxoff');
drawPublishAxis('labelFontSize', 14);

% Next plot point plots split out by texture family.
subplot(2,1,2);
handles = [];
for ti = 1:length(texFams)
  col = colors(ti,:);
  h = myerrorbar((1:length(rois))-0.1, corrs.('within_texnoise')(:,ti), 'yError', corrs.('within_textex_se')(:,ti), 'Symbol', 'o', 'Color', col, 'MarkerSize', 10); hold on;
  handles = [handles h];
  h = myerrorbar((1:length(rois))+0.1, corrs.('between_texnoise')(:,ti), 'yError', corrs.('between_textex_se')(:,ti), 'Symbol', 'o', 'Color', col+(1-col)*(.5), 'MarkerSize', 10); hold on;
  handles = [handles h];
end
ylim([-.1, .3]);
set(gca, 'XTick', 1:length(rois));
set(gca, 'XTickLabel', rois);
title('Correlation between textures and phase-scrambles (within vs between family), split by texture family');
xlim([0 length(rois)+1]);
labels=[cellfun(@(x) [x ' Within'], stimTraces.imNames, 'un', 0), cellfun(@(x) [x ' Between'], stimTraces.imNames, 'un', 0)];
legend(handles([1:2:end, 2:2:end]), labels, 'Location', 'Northwest');
legend('boxoff');
drawPublishAxis('labelFontSize', 14);

set(f, 'Position', [43.4, 16, 23.1, 21.2]);

if saveFigs
  savepdf(f, [saveLoc '/withinbetween_texnoise_correlation.pdf']);
end

%% Next run within/between correlation analysis, with different numbers of voxels.
nVoxels = [25, 50, 100, 150, 200, 350, 500, 750, 1000];
within_corr = zeros(length(rois), length(texFams), length(nVoxels));
between_corr = zeros(length(rois), length(texFams), length(nVoxels));
parfor i = 1:length(nVoxels)
  nVox = nVoxels(i);
  nCorr_mtx = ER_getCorrMtx(rois, stimTraces, [], nVox);
  nCorrs = ER_getCorrStruct(rois, stimTraces, nCorr_mtx);

  within_corr(:,:,i) = nCorrs.within_textex;
  within_se(:,:,i) = nCorrs.within_textex_se;
  between_corr(:,:,i) = nCorrs.between_textex;
  between_se(:,:,i) = nCorrs.between_textex_se;
end
% Then, plot within/between family correlation as a function of number of voxels used in the analysis.
wc = squeeze(mean(within_corr, 2)); % nROIs x length(nVoxels)
bc = squeeze(mean(between_corr, 2));

f=figure;
handles = [];
for i = 1:size(wc,1)
  col = colors(i,:);
  h1 = plot(nVoxels, wc(i,:), '.-', 'Color', colors(i,:)); hold on;
  text(nVoxels(end), wc(i,end), sprintf('%s: within', rois{i}));
  handles = [handles h1];
end
for i = 1:size(wc,1);
  col = colors(i,:);
  h2 = plot(nVoxels, bc(i,:), '.-', 'Color', col+(1-col)*(.5)); 
  text(nVoxels(end), bc(i,end), sprintf('%s: between', rois{i}))
  handles = [handles h2];
end

xlabel('Number of voxels included in analysis');
ylabel('Split-half correlation');
title('Split-half correlation as a function of number of voxels');
legend(handles, [cellfun(@(x) [x ': within'], rois, 'un', 0) cellfun(@(x) [x ': between'], rois, 'un', 0)]);
legend('boxoff');
drawPublishAxis('labelFontSize', 14);
set(f, 'Position', [23.8, 13.9, 33.5, 23.1]);

savepdf(f, [saveLoc '/splithalf_vs_nVoxels']);
saveas(f, [saveLoc '/splithalf_vs_nVoxels.png']);

%% Next run within/between correlation analysis, with different r2 cutoffs.
r2s = [0, 0.05, .1, .15, .2, .25, .3, .4, .5];
within_corr = zeros(length(rois), length(texFams), length(r2s));
between_corr = zeros(length(rois), length(texFams), length(r2s));
parfor i = 1:length(r2s)
  r2i = r2s(i);
  nCorr_mtx = ER_getCorrMtx(rois, stimTraces, r2i);
  nCorrs = ER_getCorrStruct(rois, stimTraces, nCorr_mtx);

  within_corr(:,:,i) = nCorrs.within_textex;
  within_se(:,:,i) = nCorrs.within_textex_se;
  between_corr(:,:,i) = nCorrs.between_textex;
  between_se(:,:,i) = nCorrs.between_textex_se;
end
% Then, plot within/between family correlation as a function of number of voxels used in the analysis.
wc = squeeze(mean(within_corr, 2)); % nROIs x length(nVoxels)
bc = squeeze(mean(between_corr, 2));

f=figure;
handles = [];
for i = 1:size(wc,1)
  col = colors(i,:);
  h1 = plot(r2s, wc(i,:), '.-', 'Color', colors(i,:)); hold on;
  text(r2s(end), wc(i,end), sprintf('%s: within', rois{i}));
  handles = [handles h1];
end
for i = 1:size(wc,1);
  col = colors(i,:);
  h2 = plot(r2s, bc(i,:), '.-', 'Color', col+(1-col)*(.5)); 
  text(r2s(end), bc(i,end), sprintf('%s: between', rois{i}))
  handles = [handles h2];
end

xlabel('r2 threshold to include voxels');
ylabel('Split-half correlation');
title('Split-half correlation as a function of r2 threshold');
legend(handles, [cellfun(@(x) [x ': within'], rois, 'un', 0) cellfun(@(x) [x ': between'], rois, 'un', 0)]);
legend('boxoff');
drawPublishAxis('labelFontSize', 14);
set(f, 'Position', [23.8, 13.9, 33.5, 23.1]);

savepdf(f, [saveLoc '/splithalf_vs_r2']);
saveas(f, [saveLoc '/splithalf_vs_r2.png']);

%%
keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = appendToStruc(structure, field, vector)

if isfield(structure, field)
  retval = [structure.(field) vector];
else
  retval = vector;
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to compute modulation index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mod_index, mod_index_CI] = ER_computeModulationIndex(rois, stimTraces, r2cutoff)

% Extract conditions and values
stimValues = stimTraces.stimValues; % 2x21 array
condNames = stimTraces.condNames; % Conditions: Texture Family x Sample Index

texIdx = [1 2 3 4]; % which samples are textures (textures = 1-4, noise = 5-7)
isTexture = arrayfun(@(x) any(x==texIdx), stimValues(2,:));

texFams = unique(stimValues(1,:));

mod_index = zeros(length(rois), length(texFams));
mod_index_CI = zeros(length(rois), length(texFams));

% Compute modulation index for each ROI
for ri = 1:length(rois)
  roi = rois{ri};
  
  %roiTrace = stimTraces.(sprintf('%s_traces', roi)); % 1x21 cell array
  roiAmps = stimTraces.(roi).amplitudes;
  roiR2 = stimTraces.(roi).r2;
  
  % Split out by each texture class.
  for tfI = 1:length(texFams)
    thisFam = (stimValues(1,:) == tfI);
    
    % Get the traces of this texture family + texture sample
    thisFamTex = (thisFam & isTexture);
    thisFamTexAmps = roiAmps(thisFamTex,roiR2>r2cutoff);
    
    thisFamNoise = (thisFam & ~isTexture);
    thisFamNoiseAmps = roiAmps(thisFamNoise,roiR2>r2cutoff);
    
    % Get the peak amplitude averaged across traces.
    thisTexAmp = mean(thisFamTexAmps, 1); % 1 x nVoxels
    thisNoiseAmp = mean(thisFamNoiseAmps, 1); % 1 x nVoxels
    
    %thisTexAmp = squeeze(max(mean(thisFamTexTraces, 3), [], 2)); % nVoxels x 1
    %thisNoiseAmp = squeeze(max(mean(thisFamNoiseTraces, 3), [], 2)); % nVoxels x 1
    
    % Compute modulation index for each voxel:
    this_mod_idx = (thisTexAmp - thisNoiseAmp) ./ (thisTexAmp + thisNoiseAmp);
    mod_index(ri, tfI) = mean(this_mod_idx);
    mod_index_CI(ri,tfI) = 1.96*std(this_mod_idx)/sqrt(length(this_mod_idx));   
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get full correlation matrix
% - can call either with a r2cutoff or a number of voxels.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function corr_mtx = ER_getCorrMtx(rois, stimTraces, r2cutoff, nVoxels)

cond_pairs = combvecs(1:length(stimTraces.stimNames), 1:length(stimTraces.stimNames));
corr_mtx = zeros(length(rois), length(stimTraces.stimNames), length(stimTraces.stimNames));
corr_cond = zeros(length(rois), 4);

nSmps = 11;

if ieNotDefined('r2cutoff')
  r2cutoff = .2;
end

for ri = 1:length(rois)
  traces = stimTraces.(rois{ri}).traces;
  r2s = stimTraces.(rois{ri}).r2;
  % Get voxel indices sorted by r2s.
  [sortedR2, sortIdx] = sort(r2s, 'descend');
  
  if ieNotDefined('nVoxels')
    % use r2cutoff to determine number of voxels.
    a = find(sortedR2<r2cutoff);
    if isempty(a) % all voxels exceed r2 threshold
      nVox = length(r2s);
    else
      nVox = a(1)-1;
    end
    if nVox == 1, nVox = 0; end
  else
    %
    nVox = min(length(r2s), nVoxels);
  end
  
  disppercent(-inf, sprintf('Correlation Matrix of ROI: %s', rois{ri}));
  for ci = 1:size(cond_pairs,1)
    cond = cond_pairs(ci,:);

    % first average over time. (nVoxels x nRepeats)
    t1 = squeeze(mean(traces{cond(1)},2));
    t2 = squeeze(mean(traces{cond(2)},2));

    % Draw 1000 samples of prespecified size and take correlations.
    corrs = [];
    for i = 1:100
      idx = randperm(min(size(t2,2), size(t1,2))); % randomly permute the number of observations

      % Take the average over nSmps presentations, and sort voxels by r2.
      a1 = mean(t1(sortIdx, idx(1:nSmps)),2); % take nSmps observations of condition 1
      a2 = mean(t2(sortIdx, idx(nSmps+1 : 2*nSmps)), 2); % and a different nSmps observations of condition 2

      % Take the nVox voxels with highest r2.
      c = corrcoef(a1(1:nVox), a2(1:nVox)); % and correlate them.
      corrs(i) = c(1,2);
    end
    corr_mtx(ri, cond(1), cond(2)) = mean(corrs);

    disppercent(ci / size(cond_pairs,1));
  end
  disppercent(inf);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get within/between family correlations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function corrs = ER_getCorrStruct(rois, stimTraces, corr_mtx)

% Extract conditions and values
stimValues = stimTraces.stimValues; % 2x21 array
condNames = stimTraces.condNames; % Conditions: Texture Family x Sample Index

texIdx = [1 2 3 4]; % which samples are textures (textures = 1-4, noise = 5-7)
isTexture = arrayfun(@(x) any(x==texIdx), stimValues(2,:));

texFams = unique(stimValues(1,:));

% 
corrs = {};
for ri = 1:length(rois)
  cm = squeeze(corr_mtx(ri,:,:));
  
  roi = {};
  for ti = 1:size(cm,1)
    cmI = cm(ti,:);
    thisFam = stimTraces.stimValues(1,ti);
    thisSmp = stimTraces.stimValues(2,ti);
    thisSmpType = isTexture(ti);
    
    sameFam = (stimTraces.stimValues(1,:) == thisFam);
    sameSmpType = (isTexture == thisSmpType);

    if thisSmpType == 1
      roi.(sprintf('fam%i_within_textex', thisFam)) = appendToStruc(roi,sprintf('fam%i_within_textex', thisFam), cmI(sameFam & sameSmpType));
      roi.(sprintf('fam%i_between_textex', thisFam)) = appendToStruc(roi,sprintf('fam%i_between_textex', thisFam), cmI(~sameFam & sameSmpType));
    
      roi.(sprintf('fam%i_within_texnoise', thisFam)) = appendToStruc(roi,sprintf('fam%i_within_texnoise', thisFam), cmI(sameFam & ~sameSmpType));
      roi.(sprintf('fam%i_between_texnoise', thisFam)) = appendToStruc(roi,sprintf('fam%i_between_texnoise', thisFam), cmI(~sameFam & ~sameSmpType));
    end
    
  end
  corrs.(rois{ri}) = roi;
end

corrTypes = {'within_textex', 'between_textex', 'within_texnoise', 'between_texnoise'};

for ci = 1:length(corrTypes)
  corrType = corrTypes{ci};
  corrs.(corrType) = zeros(length(rois), length(texFams));
  corrs.(sprintf('%s_se', corrType)) = zeros(length(rois), length(texFams));
  for ri = 1:length(rois)
    corrStruc = corrs.(rois{ri});
    for i = 1:length(texFams)
      cor = corrStruc.(sprintf('fam%i_%s', i, corrType));
      corrs.(corrType)(ri, i) = mean(cor);
      corrs.(sprintf('%s_se', corrType))(ri,i) = 1.96*std(cor)./length(cor);    
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get all combinations of two vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cv = combvecs(A,B)
[a,b] = meshgrid(A,B);
c = cat(2,a',b');
cv = reshape(c, [], 2);
