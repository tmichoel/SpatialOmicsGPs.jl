# Define simple test data
n = 4;
K = [1. .5 .25 .125;
    .5 1. .5 .25;
    .25 .5 1. .5;
    .125 .25 .5 1]; # kernel matrix
X = [1. 1.; 1. 1.; 1. 0; 1. 0.] # fixed effects matrix
y = [0.25; 0.5; 0.25; 0.125] # response vector

# Do the eigendecomposition of the kernel matrix and rotation of the response vector for zero covariates
λ0, U0 = eigen(K);
yr0 = U0' * y;

# Do the eigendecomposition of the kernel matrix and rotation of the response vector for one covariate
y1, K1 = project_orth_X(y, K, X[:,1]);
λ1, U1 = eigen(K1);
yr1 = U1' * y1;

# Do the eigendecomposition of the kernel matrix and rotation of the response vector for two covariates
y2, K2 = project_orth_X(y, K, X);
λ2, U2 = eigen(K2);
yr2 = U2' * y2;

@testset "FaSTLMM full-rank tests" begin 
    # set a value for δ
    δ = 0.5
    
    # Compute the MLE of the variance using FastLMM and zero covariates
    σ²0 = sigma2_mle_fullrank(δ, λ0, yr0)
    # Compute the MLE of the variance using the exact formula
    σ²0_exact = dot(y,  inv(K + δ*I) * y) / n
    # Compare the two
    @test σ²0 ≈ σ²0_exact atol=1e-6

    # Do the same for one covariate
    σ²1 = sigma2_mle_fullrank(δ, λ1, yr1)
    yhat1 = X[:,1]*inv(X[:,1]'*inv(K + δ*I)*X[:,1])*X[:,1]'*inv(K + δ*I)*y
    σ²1_exact = dot((y - yhat1), inv(K + δ*I) * (y - yhat1)) / (n - 1)
    @test σ²1 ≈ σ²1_exact atol=1e-6

    # Do the same for two covariates
    σ²2 = sigma2_mle_fullrank(δ, λ2, yr2)
    yhat2 = X*inv(X'*inv(K + δ*I)*X)*X'*inv(K + δ*I)*y
    σ²2_exact = dot((y - yhat2), inv(K + δ*I) * (y - yhat2)) / (n - 2)
    @test σ²2 ≈ σ²2_exact atol=1e-6
    
    # Compute the minus log-likelihood using FastLMM and zero covariates
    log_like0 = minus_log_like_fullrank(δ, λ0, yr0)
    # Compute the minus log-likelihood using the exact formula
    log_like0_exact = log(det(K + δ*I)) / n + log(σ²0_exact)
    # Compare the two
    @test log_like0 ≈ log_like0_exact atol=1e-6

    # Do the same for one covariate
    log_like1 = minus_log_like_fullrank(δ, λ1, yr1)
    log_like1_exact = log(det(K + δ*I)) / (n-1) + log(det(X[:,1]' * inv(K + δ*I) * X[:,1])) / (n-1) - log(dot(X[:,1],X[:,1])) / (n-1) + log(σ²1_exact)
    @test log_like1 ≈ log_like1_exact atol=1e-6

    # Do the same for two covariates
    log_like2 = minus_log_like_fullrank(δ, λ2, yr2)
    log_like2_exact = log(det(K + δ*I)) / (n-2) + log(det(X' * inv(K + δ*I) * X)) / (n-2) - log(det(X'*X)) / (n-2) + log(σ²2_exact)
    @test log_like2 ≈ log_like2_exact atol=1e-6

    # Test if the minus log-likelihood function works on a vector of inputs
    δs = [0.1, 0.5, 1.0]
    f(δ) = minus_log_like_fullrank(δ, λ0, yr0)
    log_likes = f.(δs)
    log_likes_exact = [minus_log_like_fullrank(δ, λ0, yr0) for δ in δs]
    @test log_likes ≈ log_likes_exact atol=1e-6
 end

