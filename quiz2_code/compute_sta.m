% function [ sta ] = compute_sta( stim, rho, num_timesteps )

function [ sta ] = compute_sta( stim, rho, num_timesteps )




%COMPUTE_STA Calculates the spike-triggered average for a neuron that
%            is driven by a stimulus defined in stim. The spike-
%            triggered average is computed over num_timesteps timesteps.
    sta = zeros(num_timesteps, 1);

    % This command finds the indices of all of the spikes that occur
    % after 300 ms into the recording.
    spike_times = find(rho(num_timesteps+1:end)) + num_timesteps;

    % Fill in this value. Note that you should not count spikes that occur
    % before 300 ms into the recording.

    num_spikes = size(spike_times);

    % Compute the spike-triggered average of the spikes found using the
    % find command. To do this, compute the average of all of the vectors
    % starting 300 ms (exclusive) before a spike and ending at the time of
    % the event (inclusive). Each of these vectors defines a list of
    % samples that is contained within a window of 300 ms before the each
    % spike. The average of these vectors should be completed in an
    % element-wise manner.
    %
    % Your code goes here.






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE A NEW MATRIX TO HOLD TEMPORARY VALUES

data_points = zeros(num_spikes(1), num_timesteps);

% INITIALIZE I, J

i = 1;

j = 1;

% VISIT EACH SPIKE (COLUMN)

% while j <= num_spikes(1)
while j <= 600

        % COLLECT DATA FROM EACH SPIKE
        % ... BY PICKING UP THE VALUES 149, 148, ... , 2, 1 STEPS BEFORE

        while i <= num_timesteps
        % while i <= 5

            % ADD THEM TO THE TEMPORARY MATRIX

            data_points(i,j) = stim(spike_times(j) - num_timesteps + i);
            % data_points(i,j) = stim(spike_times(j) - 5 + i);

            % PICK UP NEXT VALUE (NEXT ROW, SAME COLUMN)

            i += 1;

        endwhile


        % RESET i TO MOVE TO THE TOP OF THE MATRIX

        i = 1;

% MOVE TO NEXT SPIKE (COLUMN)

    j += 1;


% END THE WHILE LOOP

endwhile




% TAKE THE AVERAGE OF EACH ROW IN THE TEMPORARY MATRIX
% ASSIGN THAT (AVERAGE) VALUE TO THE "STA" MATRIX



% M = mean(A,dim) returns the mean values for elements along the dimension of A specified by scalar dim. For matrices, mean(A,2) is a column vector containing the mean value of each row.
% http://www.mathworks.com/help/matlab/ref/mean.html

% sta = mean(data_points, 2);

k = 1;
while k <= num_timesteps
    sta(k, 1) = mean(data_points(k, :));
    k += 1;
endwhile




% "quiz.m" WILL TAKE CARE OF PLOTTING & THE REST


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








end

