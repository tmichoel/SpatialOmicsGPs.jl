"""
    fastlmm_fullrank(y,X,K)

For a linear mixed  / Gaussian process with a full-rank kernel matrix,

```math
y \sim N\bigl(X\beta, \sigma^2(K + \delta I) \bigr)
``

where ``y`` is the response vector, ``X`` is a matrix of covariates, and ``K`` is the kernel matrix, compute the MLEs of the fixed effects weights, the variance parameter and the variance ratio using FaST-LMM with a full-rank kernel matrix.
"""
function fastlmm_fullrank(y,X,K)
    # Eigendecomposition of the kernel matrix
    EF = eigen(K);
    λ = EF.values;
    U = EF.vectors;

    # Rotate the data
    Xr = U' * X;
    yr = U' * y;

    # Compute the MLE of the variance ration δ
    δ, res = delta_mle_fullrank(λ, yr, Xr);

    # Compute the MLE of the fixed effects weights given δ
    β = beta_mle_fullrank(δ, λ, yr, Xr);

    # Compute the MLE of the variance parameter given δ and β
    σ² = sigma2_mle_fullrank(δ, λ, β, yr, Xr);

    # return the MLEs and the optimization result
    return β, σ², δ, res
end

"""
    delta_mle_fullrank(λ, yr, Xr=[])

Compute the MLE of the variance ratio δ given the eigenvalues of the kernel matrix, the rotated response vector and (optionally) the rotated covariates, by solving the non-convex optimization problem formulated in the FaST-LMM paper.
"""
function delta_mle_fullrank(λ, yr, Xr=[])
    # Compute the MLE of the variance parameter
    res = Optim.optimize(x -> minus_log_like_fullrank(softplus(x), λ, yr, Xr), [0.0], LBFGS(); autodiff = :forward);
    xmin = res.minimizer;
    return softplus(xmin), res
end

"""
    minus_log_like_fullrank(δ, λ, yr, Xr=[])

Compute the minus log-likelihood of the model, scaled by the number of samples and without constant factors, given the variance ratio δ, the eigenvalues of the kernel matrix, the rotated response vector and (optionally) the rotated covariates.
"""
function minus_log_like_fullrank(δ, λ, yr, Xr=[])
    # Compute the minus log-likelihood of the model, scaled by the number of samples and without constant factors
    #δ = softplus(x);
    β = beta_mle_fullrank(δ, λ, yr, Xr);
    σ² = sigma2_mle_fullrank(δ, λ, β, yr, Xr);
    return  0.5 * mean(log.(λ .+ δ)) .+ 0.5 * log(σ²)
end

"""
    sigma2_mle_fullrank(δ, λ, β, yr, Xr=[])

Compute the MLE of the variance parameter given the variance ratio δ, the eigenvalues of the kernel matrix, the fixed effects weights, the rotated response vector and (optionally) the rotated covariates.
"""
function sigma2_mle_fullrank(δ, λ, β, yr, Xr=[])
    # Compute the MLE of the variance parameter
    if isempty(Xr)
        return mean(yr.^2 ./ (λ .+ δ))
    else
        return mean( (yr .- Xr * β).^2 ./ (λ .+ δ) )
    end
end

"""
    beta_mle_fullrank(δ, λ, yr, Xr=[])

Compute the MLE of the fixed effects weights given the variance ratio δ, the eigenvalues of the kernel matrix, the rotated response vector and (optionally) the rotated covariates.
"""
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

"""
    softplus(x)

Compute the softplus function, i.e. log(1 + exp(x)), element-wise.
"""
function softplus(x)
    return log.(1 .+ exp.(x))
end

