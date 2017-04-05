function main()
    % Find all PCG recordings in the validation dataset
    PCG_Files = dir('../validation_dataset/*.wav');
    %% Load the trained parameter matrices for Springer's HSMM model.
    % The parameters were trained using 409 heart sounds from MIT heart
    % sound database, i.e., recordings a0001-a0409.
    load('Springer_B_matrix.mat');
    load('Springer_pi_vector.mat');
    load('Springer_total_obs_distribution.mat');

    % Load default options for classification
    springer_options   = default_Springer_HSMM_options;

    % Generate s1 and s2 heart sound segmentations for each recording
    for i=1:numel(PCG_Files)

        PCGPath = fullfile(PCG_Files(i).folder, PCG_Files(i).name)
        % Read in PCG data
        [PCG, Fs1] = audioread(PCGPath);  % load data

        % Resample data
        PCG_resampled      = resample(PCG,springer_options.audio_Fs,Fs1); % resample to springer_options.audio_Fs (1000 Hz)

        % Running runSpringerSegmentationAlgorithm.m to obtain the
        % assigned_states
        [segmentations] = runSpringerSegmentationAlgorithm(PCG_resampled, springer_options.audio_Fs, Springer_B_matrix, Springer_pi_vector, Springer_total_obs_distribution, false); % obtain the locations for S1, systole, s2 and diastole

        % Remove systole and diastole segmentations as these will not be
        % created in the real-time implementation
        segmentations(segmentations == 4) = 3;
        segmentations(segmentations == 2) = 1;
        segmentations(segmentations == 3) = 2;

        % Plot segmentation against resampled PCG data
        % plot(PCG_resampled);
        % hold on;
        % plot(segmentations);

        % Find changes in beat type classification
        changeIndx = find(abs(diff(segmentations)))+1;
        % Get type of new beat classification (1 = first heart sound (lub), 2 = 2nd
        % heart sound (dub)
        segVals = segmentations(changeIndx);

        % Save data to CSV files
        csvdata = [changeIndx, segVals];
        [~, basename, ~] = fileparts(PCGPath)
        csvpath = fullfile('../SpringerSegmentationData', strcat(basename, '_segs.csv'));
        csvwrite(csvpath, csvdata);
    end
end

