# Define simple test data
n = 3;
K = [1. .5 .25;.5 1. .5; .25 .5 1]; # kernel matrix
X = [1. 1.; 1. 1.; 1. 0] # fixed effects matrix
y = [0.25; 0.5; 0.25] # response vector

# Eigendecomposition of the kernel matrix
EF = eigen(K);

# Rotate the data
Xr = EF.vectors' * X;
yr = EF.vectors' * y;

@testset "FaSTLMM full-rank tests" begin 
    # set a value for δ
    δ = 0.5
    # Compute the MLE of the fixed effects weights using FastLMM
    β = beta_mle_fullrank(δ, EF.values, yr, Xr)
    # Compute the MLE of the fixed effects weights using the exact formula
    β_exact = inv(X'*inv(K + δ*I)*X)*X'*inv(K + δ*I)*y
    # Compare the two
    @test β ≈ β_exact atol=1e-6
    
    # Compute the MLE of the variance using FastLMM
    σ² = sigma2_mle_fullrank(δ, EF.values, β, yr, Xr)
    # Compute the MLE of the variance using the exact formula
    yhat = X*β_exact
    σ²_exact = dot((y - yhat),  inv(K + δ*I) * (y - yhat)) / n
    # Compare the two
    @test σ² ≈ σ²_exact atol=1e-6
    
    # Compute the minus log-likelihood using FastLMM
    log_like = minus_log_like_fullrank(δ; λ=EF.values, yr=yr, Xr=Xr)
    # Compute the minus log-likelihood using the exact formula
    log_like_exact = 0.5 * log(det(K + δ*I)) / n + 0.5 * log(σ²_exact)
    # Compare the two
    @test log_like ≈ log_like_exact atol=1e-6
 end

