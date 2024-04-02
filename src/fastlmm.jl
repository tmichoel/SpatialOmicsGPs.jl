function fastlmm_fullrank(y,X,K)
    # Eigendecomposition of the kernel matrix
    EF = eigen(K);
    λ = EF.values;
    U = EF.vectors;

    # Rotate the data
    Xr = U' * X;
    yr = U' * y;

    # Compute the MLE of the fixed effects weights and variance parameters
    β, σ², δ = fastlmm_mle(yr, Xr, λ);
end

function delta_mle_fullrank(λ, yr, Xr=[])
    # Compute the MLE of the variance parameter
    res = Optim.optimize(x -> minus_log_like_fullrank(softplus(x), λ, yr, Xr), [0.0], LBFGS(); autodiff = :forward);
    xmin = res.minimizer;
    return softplus(xmin), res
end

function minus_log_like_fullrank(δ, λ, yr, Xr=[])
    # Compute the minus log-likelihood of the model, scaled by the number of samples and without constant factors
    #δ = softplus(x);
    β = beta_mle_fullrank(δ, λ, yr, Xr);
    σ² = sigma2_mle_fullrank(δ, λ, β, yr, Xr);
    return  0.5 * mean(log.(λ .+ δ)) .+ 0.5 * log(σ²)
end

function sigma2_mle_fullrank(δ, λ, β, yr, Xr=[])
    # Compute the MLE of the variance parameter
    if isempty(Xr)
        return mean(yr.^2 ./ (λ .+ δ))
    else
        return mean( (yr .- Xr * β).^2 ./ (λ .+ δ) )
    end
end

function beta_mle_fullrank(δ, λ, yr, Xr=[])
    # Compute the MLE of the fixed effects weights
    if isempty(Xr)
        return 0.0
    else
        A = sum([Xr[i,:] * Xr[i,:]' ./ (λ[i] + δ) for i in eachindex(λ)]) # counterintuitive order of transpose in Xr, because a Matrix row is a (column) Vector
        b = sum([Xr[i,:] * yr[i] ./ (λ[i] + δ) for i in eachindex(λ)]) # counterintuitive order of transpose in Xr, because a Matrix row is a (column) Vector
        return A \ b
    end  
end

function softplus(x)
    return log.(1 .+ exp.(x))
end

