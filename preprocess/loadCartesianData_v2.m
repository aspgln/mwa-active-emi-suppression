function dataObj = loadCartesianData(dat_dir,dat_filename, rparams)
    % LOADCARTESIANDATA Load and preprocess Cartesian k-space data
    %
    % Inputs:
    %   twixObj - Twix object containing MRI data
    %   rparams - Configuration structure with fields:
    %             .decorrelation - Enable decorrelation (0=skip, 1=enable)
    %             .compressionNum - Number of coils for compression (0=skip, >0=compress to N coils)
    %             .loadNoise - Load pre-scan noise data (default: true)
    %             .loadACS - Load ACS reference data (default: true)
    %             .verbose - Display progress messages (default: true)
    %
    % Output:
    %   dataObj - Structure containing processed data and metadata

    % Set default parameters
    rparams = setDefaultParams(rparams);
    
    % Initialize data object
    dataObj = initializeDataObject();
    
    % Load twix object
    twixObj = mapVBVD([dat_dir,dat_filename]);

    % Load noise data
    if rparams.loadNoise
        noise = loadNoiseData(twixObj, rparams.verbose);
    else
        noise = [];
    end
    
    % Load k-space data
    kdata = loadKspaceData(twixObj, rparams.verbose);
    
    % Load ACS data
    [acs, acsFLAG] = loadACSData(twixObj, rparams.loadACS, rparams.verbose);
    
    % ====================================================================
    % PREPROCESSING STEPS - Add/remove/reorder as needed
    % ====================================================================
    
    % Step 1: Signal decorrelation
    if rparams.decorrelation ~= 0
        if rparams.verbose
            fprintf('Applying decorrelation... ');
            tic;
        end
        kdata = signal_decorrelation(kdata, noise);
        toc
        if ~isequal(acs, [0])
            acs = signal_decorrelation(acs, noise);
        end
        if rparams.verbose
            fprintf('Done. ');
            toc;
        end
    end
    
    % Step 2: Coil compression
    if rparams.compressionNum ~= 0
        if rparams.verbose
            fprintf('Applying coil compression to %d coils... ', rparams.compressionNum);
            tic;
        end
        kdata = preprocess.util.coil_selection(kdata, 'svd', rparams.compressionNum);
        if rparams.verbose
            fprintf('Done. ');
            toc;
        end
    end
    
    % % Step 3: Phase correction (example - add your own)
    % if isfield(rparams, 'phaseCorrection') && rparams.phaseCorrection
    %     if rparams.verbose
    %         fprintf('Applying phase correction... ');
    %         tic;
    %     end
    %     kdata = preprocess.phase_correction(kdata);
    %     if ~isequal(acs, [0])
    %         acs = preprocess.phase_correction(acs);
    %     end
    %     if rparams.verbose
    %         fprintf('Done. ');
    %         toc;
    %     end
    % end
    
    

    
    % % Step 6: Custom preprocessing (example - add your own)
    % if isfield(rparams, 'customPreprocess') && ~isempty(rparams.customPreprocess)
    %     if rparams.verbose
    %         fprintf('Applying custom preprocessing... ');
    %         tic;
    %     end
    %     kdata = preprocess.custom_function(kdata, rparams.customPreprocess);
    %     if rparams.verbose
    %         fprintf('Done. ');
    %         toc;
    %     end
    % end
    
    % ====================================================================
    % FINAL DATA ORGANIZATION
    % ====================================================================
    

    
    % Get dimension information from twix object
    sqzSize = twixObj{2}.image.sqzSize;
    sqzDims = twixObj{2}.image.sqzDims;
    
    % Store dimension information using twix naming
    for i = 1:length(sqzDims)
        switch sqzDims{i}
            case 'Col'
                dataObj.dim.Nfe = sqzSize(i);
            case 'Cha'
                dataObj.dim.Nch = sqzSize(i);
            case 'Lin'
                dataObj.dim.Npe = sqzSize(i);
            case 'Eco'
                dataObj.dim.Nec = sqzSize(i);
            case 'Rep'
                dataObj.dim.Nt = sqzSize(i);
        end
    end
    
    % Store original twix dimension info for reference
    dataObj.dim.sqzSize = sqzSize;
    dataObj.dim.sqzDims = sqzDims;
    
    % Calculate acceleration parameters
    % Note: Using kdata instead of undefined kdata_img
    [R_accel, firstPE] = recon.check_accel_ratio_v3b(kdata);
    dataObj.R_accel = R_accel;
    dataObj.firstPE = firstPE;
    
    % Store metadata
    dataObj.acsFLAG = acsFLAG;

    % convert to 5d
    if numel(size(kdata)) == 4
        kdata = util.make5d(kdata);
    elseif numel(size(kdata)) == 5
    else 
        error('Check Data Dimension!');
    end 

    % if numel(size(acs)) == 4
    %     acs = util.make5d(acs);
    % elseif numel(size(kdata)) == 5
    % elseif numel(size(kdata)) == 3
    % else 
    %     error('Check Data Dimension!');
    % end 

    % Store processed data (keep original dimension order from twix)
    dataObj.kdata = permute(kdata, [1 3 2 4 5 6]);% Nkx, Nky, Nch, Nec, Nt...
    dataObj.acs = permute(acs, [1 3 2 4 5 6]);
    dataObj.noise = noise;

end

% ========================================================================
% HELPER FUNCTIONS
% ========================================================================

function rparams = setDefaultParams(rparams)
    % Set default parameter values if not specified
    
    if ~isfield(rparams, 'decorrelation')
        rparams.decorrelation = 0;
    end
    
    if ~isfield(rparams, 'compressionNum')
        rparams.compressionNum = 0;
    end
    
    if ~isfield(rparams, 'loadNoise')
        rparams.loadNoise = true;
    end
    
    if ~isfield(rparams, 'loadACS')
        rparams.loadACS = true;
    end
    
    if ~isfield(rparams, 'verbose')
        rparams.verbose = true;
    end
end

function dataObj = initializeDataObject()
    % Initialize empty data object structure
    
    dataObj = struct();
    dataObj.kdata = [];
    dataObj.acs = [];
    dataObj.dim = struct();
    dataObj.R_accel = [];
    dataObj.firstPE = [];
end

function noise = loadNoiseData(twixObj, verbose)
    % Load pre-scan noise data
    
    if verbose
        fprintf('Loading noise data... ');
        tic;
    end
    
    noise = twixObj{1}.noise();
    
    if verbose
        fprintf('Done. ');
        toc;
    end
end

function kdata = loadKspaceData(twixObj, verbose)
    % Load k-space image data using twix dimension information
    
    if verbose
        fprintf('Loading kdata... ');
        tic;
    end
    
    % Load raw k-space data (keep original dimension order)
    kdata = squeeze(twixObj{2}.image());
    
    % Get dimension information from twix object
    sqzSize = twixObj{2}.image.sqzSize;
    sqzDims = twixObj{2}.image.sqzDims;
    
    if verbose
        fprintf('Twix dimensions: ');
        for i = 1:length(sqzDims)
            fprintf('%s=%d ', sqzDims{i}, sqzSize(i));
        end
        fprintf('Done. ');
        toc;
    end
end

function [acs, acsFLAG] = loadACSData(twixObj, loadACS, verbose)
    % Load ACS reference scan data if available and requested
    
    acs = [0];  % Default empty ACS
    acsFLAG = 0;
    
    if ~loadACS
        if verbose
            fprintf('ACS loading disabled by configuration.\n');
        end
        return;
    end
    
    if verbose
        tic;
    end
    
    if isfield(twixObj{2}, 'refscan')
        acs = squeeze(twixObj{2}.refscan());
        acsFLAG = 1;
        
        if verbose
            fprintf('ACS data loaded. ');
            toc;
        end
    else
        if verbose
            fprintf('No ACS data found in twix object.\n');
        end
    end
end