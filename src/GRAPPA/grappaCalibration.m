function [kernel_mcoil] = grappaCalibration(AtA, kSize, nCoil, coil, lambda, sampling)

if nargin < 6
	sampling = ones([kSize,nCoil]);
end

kernel_mcoil = zeros(kSize(1), kSize(2), nCoil, nCoil);

karea = kSize(1)*kSize(2);
idxY_arr = (karea+1)/2:karea:(karea*nCoil);
idxA_all = find(sampling);

for cc = 1:nCoil

    idxY = idxY_arr(cc);
    idxA = idxA_all;
    idxA(idxA==idxY) = [];
    
    Aty = AtA(idxA,idxY);
    AA = AtA(idxA,idxA);
    
    kernel = sampling*0;
    
    lambda = norm(AA,'fro')/size(AA,1)*lambda;
    
    kernel(idxA) = inv(AA + eye(size(AA))*lambda)*Aty;
    
    kernel_mcoil(:,:,:,cc) = kernel;

end
