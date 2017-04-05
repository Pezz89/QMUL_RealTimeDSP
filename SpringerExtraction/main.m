function main()
    % Find all PCG recordings in the validation dataset
    PCG_Files = dir('../validation_dataset/*.wav');
    % Generate s1 and s2 heart sound segmentations for each recording
    for i=1:numel(PCG_Files)
        segmentations = challenge(PCG_Files(i));
    end
    % Save to csv files as time, value pairs
end

