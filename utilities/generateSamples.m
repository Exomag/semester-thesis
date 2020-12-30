function samples = generateSamples(N, mean, var, dist)
    if isequal(dist, 'unif')
        b = mean + sqrt(3*var);
        a = 2 * mean - b;
        samples = a + (b - a) * rand(N, 1);
    elseif isequal(dist, 'norm')
        samples = mean + sqrt(var) * randn(N, 1);
    else
        error('unknown distribution')
    end
end
