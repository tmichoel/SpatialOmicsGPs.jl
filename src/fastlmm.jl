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

function beta_mle_fullrank(δ, λ, yr, Xr=[])
    # Compute the MLE of the fixed effects weights
    if isempty(Xr)
        return 0.0
    else
        A = sum([Xr[i,:]' * Xr[i,:] ./ (λ[i] + δ) for i in eachindex(λ)])
        b = sum([Xr[i,:]' * yr[i] ./ (λ[i] + δ) for i in eachindex(λ)])
        return A \ b
    end
    
end

function sigma2_mle_fullrank(δ, λ, β, yr, Xr=[])
    # Compute the MLE of the variance parameter
    if isempty(Xr)
        return mean(yr.^2 / (λ .+ δ))
    else
        return mean( (yr .- Xr * β).^2 ./ (λ .+ δ) )
    end
end

function log_like_fullrank(δ, λ, yr, Xr=[])
    # Compute the log-likelihood of the model
    β = beta_mle_fullrank(yr, Xr, λ, δ);
    σ² = sigma2_mle_fullrank(yr, Xr, λ, δ, β);
    return  - 0.5 * mean(log.(λ .+ δ))  - 0.5 * log(σ²)
end