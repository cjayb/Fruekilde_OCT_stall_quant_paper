function varargout = clean_angio_stripes_and_emph_vessels(img, vessel_mask)

if nargin < 2
    vessel_mask = [];
end

destripe_wavelet = 'db2';
destripe_keep_comps = 20;
destripe_sigma = 1;

% https://www.mathworks.com/matlabcentral/fileexchange/24409-hessian-based-frangi-vesselness-filter
cfg_opts = struct();
cfg_opts.FrangiScaleRange = [1, 3];  % The range of sigmas used, default [1 8]
cfg_opts.FrangiScaleRatio = 2;  % Step size between sigmas, default 2
cfg_opts.FrangiBetaOne = 0.5;  % Frangi correction constant, default 0.5
cfg_opts.FrangiBetaTwo = 10;  % Frangi correction constant, default 15
cfg_opts.BlackWhite = false;% Detect black ridges (default) set to true, for
                      % white ridges set to false.
cfg_opts.verbose = false;  % Show debug information, default true


destripe = RemoveStripesVertical(img, destripe_keep_comps, ...
                                 destripe_wavelet, destripe_sigma);    

med_destripe = medfilt2(destripe, [3 3]);

if ~isempty(vessel_mask)
    med_destripe = med_destripe .* vessel_mask;
end

frangi = FrangiFilter2D(med_destripe, cfg_opts);

clean_img = histeq(frangi / max(frangi(:)));

varargout{1} = clean_img;

% for plotting
if nargout > 1
    varargout{1} = destripe;
    varargout{2} = frangi;
    varargout{3} = clean_img;
else
    varargout{1} = clean_img;
end

%% LICENSE INFORMATION of FrangiFilter2D
% Copyright (c) 2009, Dirk-Jan Kroon
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
