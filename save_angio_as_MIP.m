function log_angio = save_angio_as_MIP(angio_raw, out_folder_tif, ...
                                       filename)

log_angio = double(log(squeeze(max(angio_raw, [], 1))));

img = log_angio;
img = (img - min(img(:))) / (max(img(:)) - min(img(:)));  % Normalization

imwrite(img, fullfile(out_folder_tif, filename));
