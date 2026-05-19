function [kdata_GRAPPA_res, img_GRAPPA_res] = GRAPPA_2D(kdata, kdata_dimStr, acs, acs_dimStr, varargin)

% Initialize opts with default values
opts.GRAPPA_processing_dim = '2D';
opts.GRAPPA_processing_contrast = 'scon';
opts.Raccel = [1,2];
opts.GRAPPAkernel = [5,5];
opts.Tk_lambda = 0.01;

% varargin handling (must be option/value pairs)
for k = 1:2:numel(varargin)
    if k==numel(varargin) || ~ischar(varargin{k})
        error('''varargin'' must be option/value pairs.');
    end
    if ~isfield(opts,varargin{k})
        warning('''%s'' is not a valid option.',varargin{k});
    end
    opts.(varargin{k}) = varargin{k+1};
end

% Load to local variables
GRAPPA_processing_dim = opts.GRAPPA_processing_dim;
GRAPPA_processing_contrast = opts.GRAPPA_processing_contrast;
Raccel = opts.Raccel;
GRAPPAkernel = opts.GRAPPAkernel;
Tk_lambda = opts.Tk_lambda;

% Adjust the "kdata" and "acs" dimension to (Nkx,Nky,Nch,Nz,Ncon,Na)
kdata_order = getGRAPPAReconDim(kdata_dimStr);
[~, kdata_order_rev] = sort(kdata_order,'ascend');
kdata = permute(kdata, kdata_order);

acs_order = getGRAPPAReconDim(acs_dimStr);
[~, acs_order_rev] = sort(acs_order,'ascend');
acs = permute(acs, acs_order);

[Nkx,Nky,Nch,Nz,Ncon,Nrpt] = size(kdata);
[Nkx_acs,Nky_acs,Nch_acs,Nz_acs,Ncon_acs,Nrpt_acs] = size(acs);

dim_check = checkGRAPPA_dimension(kdata,acs);
if dim_check == 0
    error('Error: kdata and acs dimension mismatch');
end

img_GRAPPA_res = zeros(Nkx,Nky,Nch,Nz,Ncon,Nrpt);
kdata_GRAPPA_res = zeros(Nkx,Nky,Nch,Nz,Ncon,Nrpt);
gmap_res = zeros(Nkx,Nky,Nch,Nz,Ncon,Nrpt);

if strcmp(GRAPPA_processing_contrast,'scon') == 1
    for a = 1:Nrpt
    for sl = 1:Nz
    for con = 1:Ncon
        ktmp = kdata(:,:,:,sl,con,a);
        acstmp = acs(:,:,:,sl,con,a);
        if a == 1 && sl == 1 && con == 1
            [patternId, patternList, maxPatternId] = generatePatternId(ktmp, GRAPPAkernel);
        end
        kres_tmp = grappaRecon(ktmp, acstmp, GRAPPAkernel, Tk_lambda, patternId, patternList, maxPatternId);
        kdata_GRAPPA_res(:,:,:,sl,con,a) = kres_tmp;
        img_GRAPPA_res(:,:,:,sl,con,a) = ifftshift(ifft2(fftshift(kres_tmp)));
    end
    end
    end
end





kdata_GRAPPA_res = permute(kdata_GRAPPA_res, kdata_order_rev);
img_GRAPPA_res = permute(img_GRAPPA_res, kdata_order_rev);
