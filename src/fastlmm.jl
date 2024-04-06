"""
    fastlmm_fullrank(y,K; covariates = [], mean = true)

For a linear mixed model / Gaussian process,

```math
y \\sim N\\bigl(X\\beta, \\sigma^2(K + \\delta I) \\bigr)
```

where ``y`` is the response vector, ``X`` is a matrix of covariates,  and ``K`` is  a full-rank kernel matrix, compute the REMLs of the variance parameter ``\\sigma^2`` and the variance ratio ``\\delta`` using FaST-LMM with a full-rank kernel matrix. Compared to the original FaST-LMM algorithm, we first project out the (optional) covariates, incl. an (optional) constant off-set (`mean=true`), from the response vector and the kernel matrix. This avoids all matrix computations in the variance parameter estimation, but means all variance estimates are restricted maximum-likelihood estimates (REMLs). Estimates for the fixed effects ``\\beta`` are not computed.
"""
function fastlmm_fullrank(y,K; covariates = [], mean = true)
    # Create covariate matrix X from the provided covariates with an intercept column if mean=true
    if !isempty(covariates) 
        X = covariates
        if mean
            X = hcat(ones(length(y)), X)
        end
    elseif mean
        X = ones(length(y))    
    end
    # if X is not empty, project it out from the response vector and the kernel matrix
    if !isempty(X)
        y, K = project_orth_covar(y, K, X)
    end
    
    # Eigendecomposition of the kernel matrix
    EF = eigen(K);
    λ = EF.values;
    U = EF.vectors;

    # Rotate the data
    yr = U' * y;

    # Compute the REML of the variance ration δ
    δ, res = delta_mle_fullrank(λ, yr);

    # Compute the REML of the variance parameter σ²
    σ² = sigma2_mle_fullrank(δ, λ, yr);

    # return the MLEs and the optimization result
    return σ², δ, res
end

"""
    delta_mle_fullrank(λ, yr)

Compute the MLE of the variance ratio δ given the eigenvalues of the kernel matrix and the rotated response vector by solving the non-convex optimization problem formulated in the FaST-LMM paper.
"""
function delta_mle_fullrank(λ, yr)
    # Compute the MLE of the variance parameter
    res = Optim.optimize(x -> minus_log_like_fullrank(softplus(x), λ, yr), [0.0], LBFGS(); autodiff = :forward);
    xmin = res.minimizer;
    return softplus(xmin[1]), res
end

"""
    minus_log_like_fullrank(δ, λ, yr)

Compute the minus log-likelihood of the model, scaled by the number of samples and without constant factors, given the variance ratio δ, the eigenvalues of the kernel matrix, the rotated response vector and (optionally) the rotated covariates.
"""
function minus_log_like_fullrank(δ, λ, yr)
    # Compute the minus log-likelihood of the model, scaled by half the number of samples and without constant factors
    σ² = sigma2_mle_fullrank(δ, λ, yr);
    return  mean(log.(λ .+ δ)) .+ log(σ²)
end

"""
    sigma2_mle_fullrank(δ, λ, yr)

Compute the MLE of the variance parameter given the variance ratio δ, the eigenvalues of the kernel matrix, and the rotated response vector.
"""
function sigma2_mle_fullrank(δ, λ, yr)
    # Compute the MLE of the variance parameter
    return mean(yr.^2 ./ (λ .+ δ))
end

# """
#     beta_mle_fullrank(δ, λ, yr, Xr=[])

# Compute the MLE of the fixed effects weights given the variance ratio δ, the eigenvalues of the kernel matrix, the rotated response vector and (optionally) the rotated covariates.
# """
# function beta_mle_fullrank(δ, λ, yr, Xr=[])
#     # Compute the MLE of the fixed effects weights
#     if isempty(Xr)
#         return 0.0 # no covariates
#     elseif size(Xr,2) == 1
#         return sum(Xr .* yr ./ (λ .+ δ)) / sum(Xr.^2 ./ (λ .+ δ)) 
#     else
#         A = Xr' * diagm(1 ./ (λ .+ δ)) * Xr
#         b = Xr' * (yr ./ (λ .+ δ))
#         return A \ b
#     end  
# end

"""
    project_orth_covar(y, K, X)

Project out the covariates from the response vector and the kernel matrix.
"""
function project_orth_covar(y, K, X)
    # SVD of the covariate matrix
    U, S,  = svd(X, full=true);
    # Get the vectors that span the orthogonal complement of the column space of X
    U2 = U[:,length(S)+1:end];
    # Project out the covariates from the response vector and the kernel matrix
    y2 = U2' * y;
    K2 = U2' * K * U2;
    # make sure K2 is symmetric
    K2 = (K2 + K2') / 2
    return y2, K2
end

"""
    softplus(x)

Compute the softplus function, i.e. log(1 + exp(x)), element-wise.
"""
function softplus(x)
    return log.(1 .+ exp.(x))
end

