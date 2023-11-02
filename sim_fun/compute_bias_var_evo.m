%% Function to compute key statistics
function [x_est_all_evo, bias_all_evo, var_all_evo, std_all_evo] = compute_bias_var_evo(x_out_all_evo, points)
    x_est_all_evo = zeros(size(x_out_all_evo,2),size(x_out_all_evo,1));
    bias_all_evo = zeros(size(x_est_all_evo));
    var_all_evo = zeros(size(x_est_all_evo));
    std_all_evo = zeros(size(x_est_all_evo));

    for epoch = 1:size(x_out_all_evo, 2)
        x_out_all = x_out_all_evo(:, epoch, :);

        var_all = circ_var(x_out_all, [],[],3);
        std_all = circ_std(x_out_all, [],[],3);
        x_est_all = circ_mean(x_out_all,[],3);
        bias_all = x_est_all - points';
        bias_all(bias_all > pi) = bias_all(bias_all > pi) - 2*pi;
        bias_all(bias_all < -pi) = bias_all(bias_all < -pi) + 2*pi;

        x_est_all_evo(epoch, :) = x_est_all;
        bias_all_evo(epoch, :) = bias_all;
        var_all_evo(epoch, :) = var_all;
        std_all_evo(epoch, :) = std_all;
    end
end