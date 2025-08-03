% Script for type3 of the 2nd level delay and dispersion analysis with updated classification rule

sub = '44';
typ = '3';
n_clusters = 9;

fully_concordant = [1];
partially_concordant = [2];
partially_discordant = [8,9];
fully_discordant = [3,4,5,6,7];

path = sprintf('/home/yosra.jamiliyan1/new_data/ICE0%s/custom_basis_fun/2nd_del_dis', sub);
basis_func = load('/home/yosra.jamiliyan1/fsl/data/feat5/default_flobs.flobs/hrfbasisfns.txt');
time = linspace(0, 30, size(basis_func, 1));

groups = struct('FC', [], 'PC', [], 'PD', [], 'FD', []);
summary_all = [];

for i = 1:n_clusters
    pe7 = niftiread(fullfile(path, sprintf('pe7_cluster%d.nii.gz', i)));
    pe8 = niftiread(fullfile(path, sprintf('pe8_cluster%d.nii.gz', i)));
    pe9 = niftiread(fullfile(path, sprintf('pe9_cluster%d.nii.gz', i)));
    mask = pe7 ~= 0;
    n_voxels = nnz(mask);
    HRFs = zeros(n_voxels, size(basis_func, 1));

    v = 1;
    for x = 1:size(pe7, 1)
        for y = 1:size(pe7, 2)
            for z = 1:size(pe7, 3)
                if mask(x, y, z)
                    beta = [pe7(x, y, z), pe8(x, y, z), pe9(x, y, z)];
                    HRFs(v, :) = beta * basis_func';
                    v = v + 1;
                end
            end
        end
    end

    mean_hrf = mean(HRFs, 1);

    early_mask = (time >= 1 & time <= 9);
    late_mask = (time > 9);

    hrf_early = mean_hrf(early_mask);
    hrf_late = mean_hrf(late_mask);
    time_early = time(early_mask);
    time_late = time(late_mask);

    [pos_peaks, pos_idx] = findpeaks(hrf_early);
    [neg_peaks_early, neg_idx_early] = findpeaks(-hrf_early); neg_peaks_early = -neg_peaks_early;

    all_peaks_early = [pos_peaks, neg_peaks_early];
    all_idx_early = [pos_idx, neg_idx_early];

    if ~isempty(all_peaks_early)
        [~, early_max_idx] = max(abs(all_peaks_early));
        early_peak = all_peaks_early(early_max_idx);
        early_peak_time = time_early(all_idx_early(early_max_idx));

        [neg_peaks_late, neg_idx_late] = findpeaks(-hrf_late); neg_peaks_late = -neg_peaks_late;
        if ~isempty(neg_peaks_late)
            [neg_peak, idx_neg] = min(neg_peaks_late);
            neg_peak_time = time_late(neg_idx_late(idx_neg));
        else
            neg_peak = 0;
            neg_peak_time = NaN;
        end

        if abs(neg_peak) > 3 * abs(early_peak)
            peak_val = neg_peak;
            peak_time = neg_peak_time;
            peak_type = "negative";
            color = [0.8 0 0];
        else
            peak_val = early_peak;
            peak_time = early_peak_time;
            if peak_val > 0
                peak_type = "positive";
                color = [0 0.6 0];
            else
                peak_type = "negative";
                color = [0.8 0 0];
            end
        end
    else
        [peak_val, idx_all] = max(abs(mean_hrf));
        peak_time = time(idx_all);
        peak_val = mean_hrf(idx_all);
        if peak_val > 0
            peak_type = "positive";
            color = [0 0.6 0];
        else
            peak_type = "negative";
            color = [0.8 0 0];
        end
    end

    [~, peaks] = max(abs(HRFs), [], 2);
    peak_latencies = time(peaks);
    peak_latency_sd = std(peak_latencies);

    % FWHM (0-15s window)
    search_fwhm_mask = (time >= 0 & time <= 15);
    time_fwhm = time(search_fwhm_mask);
    hrf_fwhm = mean_hrf(search_fwhm_mask);
    half_max = peak_val / 2;

    if strcmp(peak_type, "positive")
        cross_idx = find(hrf_fwhm >= half_max);
    elseif strcmp(peak_type, "negative")
        cross_idx = find(hrf_fwhm <= half_max);
    else
        cross_idx = [];
    end

    [~, pk_idx] = min(abs(time_fwhm - peak_time));
    left = cross_idx(cross_idx < pk_idx);
    right = cross_idx(cross_idx > pk_idx);

    if ~isempty(left) && ~isempty(right)
        fwhm_peak_only = time_fwhm(right(end)) - time_fwhm(left(1));
    else
        fwhm_peak_only = NaN;
    end

    voxel_fwhms = zeros(n_voxels, 1);
    for v = 1:n_voxels
        h = HRFs(v, :);
        [maxval, ~] = max(h);
        [minval, ~] = min(h);
        if abs(maxval) >= abs(minval)
            pval = maxval; h_half = pval / 2; idx = find(h >= h_half);
        else
            pval = minval; h_half = pval / 2; idx = find(h <= h_half);
        end
        if numel(idx) >= 2
            voxel_fwhms(v) = time(idx(end)) - time(idx(1));
        else
            voxel_fwhms(v) = NaN;
        end
    end
    fwhm_sd = std(voxel_fwhms, 'omitnan');

    this_metrics = [i, peak_val, peak_time, double(peak_val > 0), peak_latency_sd, fwhm_peak_only, fwhm_sd];

    if ismember(i, fully_concordant)
        groups.FC = [groups.FC; this_metrics];
    elseif ismember(i, partially_concordant)
        groups.PC = [groups.PC; this_metrics];
    elseif ismember(i, partially_discordant)
        groups.PD = [groups.PD; this_metrics];
    elseif ismember(i, fully_discordant)
        groups.FD = [groups.FD; this_metrics];
    end

    summary_all = [summary_all; this_metrics];

    % Plot
    figure;
    plot(time, mean_hrf, 'k', 'LineWidth', 2); hold on;
    plot(peak_time, peak_val, 'o', 'MarkerSize', 6, 'MarkerFaceColor', color, 'MarkerEdgeColor', 'k');
    if ~isnan(fwhm_peak_only) && ~isempty(left) && ~isempty(right)
        plot([time_fwhm(left(1)), time_fwhm(right(end))], [half_max, half_max], 'b--', 'LineWidth', 1.5);
        plot([time_fwhm(left(1)), time_fwhm(left(1))], [0, half_max], 'b:', 'LineWidth', 1.2);
        plot([time_fwhm(right(end)), time_fwhm(right(end))], [0, half_max], 'b:', 'LineWidth', 1.2);
    end
    xlabel('Time (s)'); ylabel('BOLD Response');
    title(sprintf('Cluster %d: Averaged HRF (%s)', i, peak_type));
    legend({'Mean HRF','Peak','FWHM'}, 'Location','northeastoutside');
    drawnow;
end

output_groups = {'FC', 'PC', 'PD', 'FD'};
for g = 1:numel(output_groups)
    group = output_groups{g};
    data = groups.(group);
    if isempty(data), continue; end
    filename = sprintf('%s_int_del_dis_sub%s_type%s.txt', group, sub, typ);
    fid = fopen(filename, 'w');
    for r = 1:size(data, 1)
        fprintf(fid, '%.4f %.4f %.4f %.4f %.4f\n', data(r,2), data(r,3), data(r,5), data(r,6), data(r,7));
    end
    fclose(fid);
end

T = array2table(summary_all, 'VariableNames', {'Cluster','PeakValue','PeakTime','IsPosDominant','SD_PeakLatency','MeanFWHM','SDFWHM'});
disp(T);
