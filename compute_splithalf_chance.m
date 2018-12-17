function mean_chance_corr = compute_splithalf_chance(stimTraces)

%%
nRepeats = 10000;
disp(nRepeats);
rois = {'V1', 'V2', 'V3', 'V4'};
nObservations = [2,4,6,8,11];

chance_corrs = zeros(length(nObservations), length(rois), nRepeats);
for oi = 1:length(nObservations)
  kObs = nObservations(oi);
  for ri = 1:length(rois)
    roi_traces = cellfun(@(x) squeeze(mean(x, 2)), stimTraces.(sprintf('%s_traces', rois{ri})), 'un', 0);

    % first, stack traces from all conditions.
    all_traces = cat(2,roi_traces{:});

    for i = 1:nRepeats
      idx = randperm(size(all_traces,2)); % randomly permute the number of observations
      I_idx = idx(1:kObs); J_idx = idx(kObs+1: 2*kObs);
      % Randomly sample k observations from each conditions and correlate
      smpI = mean(all_traces(:,I_idx),2); 
      smpJ = mean(all_traces(:,J_idx),2);

      c = corrcoef(smpI, smpJ);
      corrs(i) = c(1,2);
    end
    chance_corrs(oi,ri,:) = corrs;
    
  end
end

mean_chance_corr = mean(chance_corrs, 3);



% %%
% chance_sort  = sort(chance_corrs, 3, 'descend');
% ci = 0.5*chance_sort(:,:,ceil(.025*nRepeats)) - chance_sort(:,:,floor(.975*nRepeats));
% 
% colors = brewermap(5, 'YlGnBu');
% colors = colors(2:end, :);
% 
% figure;
% for i = 1:4
%   %plot(nObservations, mean_chance_corr(:,i), '.:', 'MarkerSize', 15, 'Color', colors(i,:)); hold on;
%   myerrorbar(nObservations, mean_chance_corr(:,i), 'yError', ci(:,i), 'Symbol=:o', 'Color', colors(i,:));
% end
