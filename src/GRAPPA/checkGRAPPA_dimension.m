function dim_check = checkGRAPPA_dimension(kdata,acs)

[Nkx,Nky,Nch,Nz,Ncon,Nrpt] = size(kdata);
[Nkx_acs,Nky_acs,Nch_acs,Nz_acs,Ncon_acs,Nrpt_acs] = size(acs);

dim_check = 1;
if Nkx ~= Nkx_acs
    dim_check = 0;
    fprintf('Dimension kx of data and acs not matched\n');
end
if Nch ~= Nch_acs
    dim_check = 0;
    fprintf('Dimension ch of data and acs not matched\n');
end
if Nz ~= Nz_acs
    dim_check = 0;
    fprintf('Dimension z of data and acs not matched\n');
end
if Ncon ~= Ncon_acs
    dim_check = 0;
    fprintf('Dimension con of data and acs not matched\n');
end
if Nrpt ~= Nrpt_acs
    dim_check = 0;
    fprintf('Dimension rpt of data and acs not matched\n');
end

