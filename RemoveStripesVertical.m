function [nima] = RemoveStripesVertical(ima, decNum, wname, sigma)
% https://www.osapublishing.org/oe/fulltext.cfm?uri=oe-17-10-8567&id=179485

for ii = 1:decNum
    [ima, Ch{ii}, Cv{ii}, Cd{ii}] = dwt2(ima, wname);
end

for ii = 1:decNum
    % FFT
    fCv = fftshift(fft(Cv{ii}));
    [my, mx] = size(fCv);
    
    damp = 1 - exp(-(-floor(my/2):-floor(my/2)+my-1).^2 / (2 * sigma^2) );
    fCv = fCv .* repmat(damp', 1, mx);
    
    % iFFT
    Cv{ii} = ifft(ifftshift(fCv));
end

% WL
nima = ima;
for ii = decNum:-1:1
    nima = nima(1:size(Ch{ii}, 1), 1:size(Ch{ii}, 2));
    nima = idwt2(nima, Ch{ii}, Cv{ii}, Cd{ii}, wname);
end

