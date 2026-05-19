function [im, imCC] = perform_grappa(kdata, acs)
% Performs GRAPPA reconstruction with automatic ACS padding

[Nkx, Nky, Nch] = size(kdata);

% Zero-pad ACS if dimensions don't match kdata
if size(acs, 1) ~= Nkx
    pad_size = (Nkx - size(acs, 1)) / 2;
    acs = padarray(acs, [pad_size, 0, 0], 0, 'both');
end
        
% Grappa Recon 
[~, im] = GRAPPA_2D(kdata, 'xyc',acs, 'xyc');

% Adaptive Coil Combination  
imCC = coilCombine(im,'xyc', 'ACC', 'accel',2);
imCC = imCC/abs(max(imCC(:)));

end