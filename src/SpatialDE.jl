
"""
    spatialde_param_inference(y,x::Union{RowVecs, ColVecs},n::Int; covariates = [], mean = true, lambda_tol = 1e-3)

Perform parameter inference for the SpatialDE model for all variables (columns) in `y` observed at locations `x` with a squared exponential kernel and a grid of `n` length scales.  The input `x` should be of type `RowVecs` or `ColVecs` from the `KernelFunctions` package.  The function returns the optimal variance parameters and length scale indices for all variables, as well as the grid of length scales used. 
"""
function spatialde_param_inference(y,x::Union{RowVecs, ColVecs},n::Int; covariates = [], mean = true, lambda_tol = 1e-3)
    # Compute all variance parameter REMLs for all variables across a grid of length scales
    lσ²s, lδs, minus_loglikes, ls = spatialde_varpars_all(y, x, n; covariates = covariates, mean = mean, lambda_tol = lambda_tol)
    # Compute the optimal variance parameters and length scale indices for all variables
    σ²s, δs, lsid = spatialde_varpars_opt(lσ²s, lδs, minus_loglikes)
    # Return the results
    return σ²s, δs, lsid, ls
end

"""
    spatialde_varpars_all(y,x::Union{RowVecs, ColVecs},n::Int; covariates = [], mean = true, lambda_tol = 1e-3)

Compute the spatial variance parameters for all variables (columns) in `y` observed at locations `x` using FaST-LMM with a squared exponential kernel and a grid of length scales.  The input `x` should be of type `RowVecs` or `ColVecs` from the `KernelFunctions` package.
"""
function spatialde_varpars_all(y,x::Union{RowVecs, ColVecs},n::Int; covariates = [], mean = true, lambda_tol = 1e-3)
    # Compute the grid of length scales
    ls = lengthscales(x, n)
    # Initialize the output arrays
    lσ²s = zeros(size(y,2),length(ls))
    lδs = zeros(size(y,2),length(ls))
    minus_loglikes = zeros(size(y,2),length(ls))
    # Loop over the length scales
    for i in eachindex(ls)
        # Create the kernel matrix
        K = kernelmatrix(with_lengthscale(SqExponentialKernel(), ls[i]),x)
        # Run FaST-LMM
        lσ²s[:,i], lδs[:,i], minus_loglikes[:,i] = fastlmm_fullrank(y, K; covariates=covariates, mean=mean, lambda_tol = lambda_tol)
    end
    # Return the results
    return lσ²s, lδs, minus_loglikes, lengthscales
end

"""
    spatialde_varpars_opt(lσ²s, lδs, minus_loglikes)

Compute the optimal spatial variance parameters and length scale index for all variables by finding for each row the minimum of the negative log-likelihood in `minus_loglikes` and returning the corresponding variance parameters and length scale index.  The input `lσ²s` and `lδs` should be matrices with the same number of rows as `minus_loglikes` and the same number of columns as the length scales used in the optimization.
"""
function spatialde_varpars_opt(lσ²s, lδs, minus_loglikes)
    # Initialize the output arrays
    n = size(lσ²s,1)
    σ²s = zeros(n)
    δs = zeros(n)
    lsid = zeros(Int64,n)
    
    # Find the best length scale and corresponding variance parameters for each variable
    for i in eachindex(σ²s)
        lsid[i] = argmin(minus_loglikes[i,:])
        σ²s[i] = lσ²s[i,lsid[i]]
        δs[i] = lδs[i,lsid[i]]
    end

    # Return the results
    return σ²s, δs, lsid
end

"""
    lengthscales(X, n=10)

For a vector-of-vectors `X`, return a grid of length scales with `n` grid points, logarithmically spaced.  The boundaries of the grid are by default set as the shortest
observed distance between pairs of vectors in `X` divided by 2, and the longest observed distance`between pairs of vectors in `X` multiplied by 2.

The input `X` should be of type `RowVecs` or `ColVecs` from the `KernelFunctions` package.
"""
function lengthscales(x::Union{RowVecs, ColVecs}, n::Int=10)

    # Compute the pairwise distances
    D = pairwise(Euclidean(), x)

    # Compute the shortest and longest distances
    minD = minimum(D[D .> 0])
    maxD = maximum(D)

    # Return the grid
    return exp10.(range(log10(minD/2), stop=log10(2*maxD), length=n))
end