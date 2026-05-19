function resres = grappaRecon_sfs(kdata,acs,kernelSize,lambda, patternId, patternList, maxPatternId)

% Original from: Michael Lustig 2008
% Adapted and modified by: Shu-Fu Shih


resres = kdata*0;
AtA = acs2CaliMat(acs, kernelSize); % build coil calibrating matrix

[Nx,Ny,Nch] = size(kdata);
kData_tmp = zeroPadding(kdata,[Nx+kernelSize(1)-1, Ny+kernelSize(2)-1,Nch]);

kernel_mat = zeros(kernelSize(1), kernelSize(2), Nch, Nch, maxPatternId);
for p = 1:maxPatternId
    kernel_mat(:,:,:,:,p) = grappaCalibration(AtA,kernelSize,Nch,1,lambda,patternList(:,:,:,p));
end

kk = reshape(kernel_mat, [kernelSize(1)*kernelSize(2)*Nch,Nch,maxPatternId]);

for x = 1:Nx
    for y = 1:Ny
            
        pnum = patternId(x,y);
        tmp = kData_tmp(x:x+kernelSize(1)-1,y:y+kernelSize(2)-1,:);
        itmp = repmat(tmp(:),[1,Nch]);

        ktmp = kk(:,:,pnum);
        resres(x,y,:) = sum(ktmp.*itmp, 1);

    end
end

end