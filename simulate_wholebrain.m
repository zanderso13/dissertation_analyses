% Parameters
num_individuals = 200; % number of individuals
num_voxels = 250000; % number of voxels per individual
num_timepoints = 800; % number of time points in the BOLD signal
tr = 1; % repetition time (seconds)

% Generate and save the brain dataset
simulate_brain_dataset(num_individuals, num_voxels, num_timepoints, tr);

function simulate_brain_dataset(num_individuals, num_voxels, num_timepoints, tr)
  % Function to simulate a brain dataset with BOLD signals
  % num_individuals: number of individuals in the dataset
  % num_voxels: number of voxels per individual
  % num_timepoints: number of time points in the BOLD signal
  % tr: repetition time (time between successive MRI scans in seconds)

  % Define HRF parameters (these are typical parameters and can be adjusted)
  hrf_duration = 30; % duration of the HRF in seconds
  hrf_peak = 6; % time to peak of the HRF in seconds
  hrf_under = 16; % time to undershoot of the HRF in seconds
  hrf_dispersion = 1; % dispersion of the HRF

  % Generate the HRF using a gamma function
  t = 0:tr:hrf_duration;
  hrf = gampdf(t, hrf_peak / hrf_dispersion, hrf_dispersion) - 0.35 * gampdf(t, hrf_under / hrf_dispersion, hrf_dispersion);

  % Pre-allocate space for the dataset
  brain_data = zeros(num_individuals, num_voxels, num_timepoints);

  % Loop through each individual
  for ind = 1:num_individuals
    % Loop through each voxel
    for vox = 1:num_voxels
      % Generate random neural activity
      neural_activity = randn(1, num_timepoints);

      % Convolve neural activity with HRF to get the BOLD signal
      bold_signal = conv(neural_activity, hrf);
      bold_signal = bold_signal(1:num_timepoints); % Truncate to original length

      % Normalize the BOLD signal to have zero mean and unit variance
      bold_signal = (bold_signal - mean(bold_signal)) / std(bold_signal);

      % Store the BOLD signal in the dataset
      brain_data(ind, vox, :) = bold_signal;
    end
  end

  % Save the dataset to a file
  save('brain_data.mat', 'brain_data');
end
