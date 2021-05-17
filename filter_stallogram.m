function filt_stalls = filter_stallogram(bin_stalls, min_stall_len, bad_frames)
% filt_stalls = filter_stallogram(bin_stalls, min_stall_len)
%
% Remove stalls shorter then min_stall_len from stallogram
%
% cjb, 27 Feb 2021
% cjb, 1 May 2021: add bad_frames

n_frames = size(bin_stalls, 2);

if nargin == 3
    bin_stalls(:, bad_frames) = 0;
end

filt_stalls = bin_stalls;
[lens, vals] = runLengthEncodeRows(filt_stalls);
for ii = 1:size(filt_stalls, 1)

    runlen = lens{ii};
    values = vals{ii};

    stalls = find(values);
    stall_durations = runlen(stalls);  % length in frames of each stall
    % filter out shortest stalls
    if ~isempty(stall_durations)
        nostall_inds = find(stall_durations < min_stall_len);
        for ns = nostall_inds
            first_idx = sum(runlen(1:stalls(ns)));
            % deal with edge
            last_idx = min(first_idx + min_stall_len - 1, n_frames);
            filt_stalls(ii, first_idx:last_idx) = 0;
        end
    end
end
