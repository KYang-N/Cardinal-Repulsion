function [x_out_all, x_sens_all] = run_inference_trials_evo(x_in, con, prior, out_epochs, T, N_epoch, x, m, memnoise)
    % Runs inference trials in MATLAB.
    % Args:
    %     x_in: input stimulus orientation, scalar in radians.
    %     con: conditional probability matrix p(m|x).
    %     out_epochs: Indices of epochs to store results for.
    %     T: Number of trials.
    %     N_epoch: Number of iterations for the recursive process.
    %     x: Stimulus space.
    %     m: Internal representation in stimulus space.
    %     memnoise: Gaussian standard deviation in radians for memory noise.
    % Returns:
    %     x_out_all: Distributions of x for all epochs (not including the initial x_in) and trials.
    %     x_sens_all: Distributions of \hat{x} for all epochs and trials.

    x_out_all = zeros(length(out_epochs), T);
    x_sens_all = zeros(length(out_epochs), T);

    for trial = 1:T
        out_epoch_idx = 1;
        x_curr = x_in;

        for epoch = 1:N_epoch
            % Find particle index in table
            [~, particle] = min(abs(x - x_curr));
            prob = con(particle, :);
            m_curr = randsample(length(m), 1, true, prob);

            % Likelihood * Prior
            posterior = con(:, m_curr)' .* prior;
            x_sens = bayes_mean(x, posterior);

            % Add memory noise
            mem_noise = normrnd(0, memnoise, 1, 1);

            % Feed x_{i+1} back as input
            x_curr = mod(x_sens + mem_noise, 2*pi);

            if ismember(epoch, out_epochs)
                x_sens_all(out_epoch_idx, trial) = x_sens;
                x_out_all(out_epoch_idx, trial) = x_curr;
                out_epoch_idx = out_epoch_idx + 1;
            end
        end
    end
end

