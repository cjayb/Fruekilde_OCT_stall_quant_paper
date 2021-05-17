function varargout = save_angios_as_MIPData(angio_raw, out_folder_data, ...
    out_folder_tif, filename)

angio = double(log(squeeze(max(angio_raw, [], 1))));

img = angio;
img = (img - min(img(:))) / (max(img(:)) - min(img(:)));  % Normalization

% TRANSPOSE to make tiff conform to ML default imagesc-view
% FLIP to make tiff conform to ML default imagesc-view
% img = flip(img', 1);  % WHY the FLIP?

imwrite(img, [fullfile(out_folder_tif, filename), '_MIPz.tif']);

outname = fullfile(out_folder_data, filename);
save(outname, 'angio');        

varargout{1} = outname;