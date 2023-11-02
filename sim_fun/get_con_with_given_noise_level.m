%% Function to get conditional probability matrix
function con = get_con_with_given_noise_level(kappa, fmap_x, fmap_m, x, m)
    con = zeros(length(x), length(m));
    
    for i = 1:length(x)
        for j = 1:length(m)
            con(i, j) = exp(kappa * (cos(fmap_m(j) - fmap_x(i)) - 1));
        end
    end
    
    con = con ./ sum(con, 2);
end
