"""
    fastlmm_fullrank(y,K; covariates = [], mean = true)

For a linear mixed model / Gaussian process,

```math
y \\sim N\\bigl(X\\beta, \\sigma^2(K + \\delta I) \\bigr)
```

where ``y`` is the response vector, ``X`` is a matrix of covariates,  and ``K`` is  a full-rank kernel matrix, compute the restricted maximum-likelihood estimates (REMLs) of the variance parameter ``\\sigma^2`` and the variance ratio ``\\delta`` using FaST-LMM with a full-rank kernel matrix. Compared to the original FaST-LMM algorithm, we first project out the (optional) covariates, incl. an (optional) constant off-set (`mean=true`), from the response vector and the kernel matrix. This avoids all matrix computations in the variance parameter estimation. Estimates for the fixed effects ``\\beta`` are not computed.
"""
function fastlmm_fullrank(y,K; covariates = [], mean = true, lambda_tol = 1e-3)
    # Create covariate matrix X from the provided covariates with an intercept column if mean=true
    X = create_covariate_matrix(covariates; mean = mean, n = size(y,1))
    
    # if X is not empty, project it out from the response vector and the kernel matrix
    if !isempty(X)
        y, K = project_orth_covar(y, K, X)
    end
    
    # Eigendecomposition of the kernel matrix
    EF = eigen(K);
    λ = EF.values;
    U = EF.vectors;

    # Due to numerical issues, we set eigenvalues smaller than lambda_tol to zero
    λ[λ .< lambda_tol] .= 0.0

    # Rotate the data
    yr = U' * y;

    if size(yr,2) == 1
        # Compute the REML of the variance ratio δ
        δ, res = delta_mle_fullrank(λ, yr);
        # Compute the REML of the variance parameter σ²
        σ² = sigma2_mle_fullrank(δ, λ, yr);
        # return the MLEs and the full optimization result
        return σ², δ, res
    else 
        # Compute and return the REMLs of the variance ratio δ and variance parameter σ² for each column of yr
        δs = zeros(size(yr,2));
        σ²s = zeros(size(yr,2));
        loglikes = zeros(size(yr,2));
        for i in eachindex(axes(yr)[2])
            δs[i], res = delta_mle_fullrank(λ, yr[:,i]);
            σ²s[i] = sigma2_mle_fullrank(δs[i], λ, yr[:,i]);
            loglikes[i] = minimum(res)
        end
        # return the MLEs and the final objective values (minus log-likelihoods)
        return σ²s, δs, loglikes
    end
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
    return  mean(log.(abs.(λ .+ δ))) .+ log(σ²)
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
    beta_mle_fullrank_lazy(y, K, X, σ², δ)

Lazy implementation of the MLE of the fixed effects weights given the response vector `y`, the kernel matrix `K`, the covariates `X`, the variance parameter `σ²` and the variance ratio `δ`. This function does not use spectral factorization, and should only be applied once, using the MLEs of the variance parameters.
"""
function beta_mle_fullrank_lazy(y, K, X, σ², δ; mean = true)
    # Create the covariance matrix
    if isempty(K)
        Σ = σ² * I
    else
        Σ = σ² * (K + δ*I)
    end
    # Create the covariate matrix
    X = create_covariate_matrix(X; mean = mean, n = size(y,1))
    # Compute the MLE of the fixed effects weights
    if isempty(X)
        return 0.0 # no covariates
    elseif size(X,2) == 1
        return dot(X, Σ \ y) / dot(X, Σ \ X)
    else
        A = X' * (Σ \ X)
        b = X' * (Σ \ y) 
        return A \ b
    end  
end


"""
    create_covariate_matrix(X; mean = true, n = 1)

Create the covariate matrix from the provided covariates with an intercept column if `mean=true`.
"""
function create_covariate_matrix(X; mean = true, n = 1)
    if !isempty(X)
        n = size(X,1)
        if mean
            X = hcat(ones(n), X)
        end
    elseif mean
        X = ones(n)
    else
        X = []  
    end
    return X
end

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
    if !isempty(K)
        K2 = U2' * K * U2;
        # make sure K2 is symmetric
        K2 = (K2 + K2') / 2
    else
        K2 = []
    end
    return y2, K2
end

"""
    softplus(x)

Compute the softplus function, i.e. log(1 + exp(x)), element-wise.
"""
function softplus(x)
    return log.(1 .+ exp.(x))
end

