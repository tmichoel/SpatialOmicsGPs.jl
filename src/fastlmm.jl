function fastlmm(y,X,K)
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