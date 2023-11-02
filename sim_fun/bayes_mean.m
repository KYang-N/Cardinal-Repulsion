function bayes_mean = bayes_mean(x, posterior)
    % Computes Bayesian estimate, minimizing L2-loss.
    num_sin = sum(sin(x) .* posterior);
    num_cos = sum(cos(x) .* posterior);
    bayes_mean = atan2(num_sin, num_cos);
    bayes_mean = bayes_mean - floor(bayes_mean / (2*pi)) * 2*pi;
end