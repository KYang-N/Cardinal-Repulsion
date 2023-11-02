%% Define cumulative function
function cumulative = cumulative(j, length, prior)
    i = round(j/length*(length));
    cumulative = sum(prior(1:i)) / sum(prior);
end