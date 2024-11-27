function momentConfig(X,U)
    Z=[vcat(X[i],U[i]) for i in eachindex(X)]
    meanZ = mean(Z)
    covZ = cov(Z)
    return [meanZ[1],meanZ[2],meanZ[3],meanZ[4],covZ[1,1],covZ[1,2],covZ[1,3],covZ[1,4],covZ[2,2],covZ[2,3],covZ[2,4],covZ[3,3],covZ[3,4],covZ[4,4]]
end

function AccMomentConfig(X,U,A)
    Z=[vcat(X[i],U[i],A[i]) for i in eachindex(X)]
    meanZ = mean(Z)
    covZ = cov(Z)
    return [meanZ[1],meanZ[2],meanZ[3],meanZ[4],meanZ[5],meanZ[6],covZ[1,1],covZ[1,2],covZ[1,3],covZ[1,4],covZ[1,5],covZ[1,6],covZ[2,2],covZ[2,3],covZ[2,4],covZ[2,5],covZ[2,6],
    covZ[3,3],covZ[3,4],covZ[3,5],covZ[3,6],covZ[4,4],covZ[4,5],covZ[4,6],covZ[5,5],covZ[5,6],covZ[6,6]]
end

function zScores(X)
    meanX = mean(X)
    varX = var(X)
    return [(x-meanX)/varX for x in X]
end

function zScoreFeatures(X,U)
    zX = zScores(X)
    zU =zScores(U)
    return momentConfig(zX,zU)
end

function zScoreFeatures(X,U,A)
    zX = zScores(X)
    zU = zScores(U)
    zA = zScores(A)
    return AccMomentConfig(zX,zU,zA)
end

function DiffMap(D,ϵ;maxoutdim=2,maxiter=1000)
    Dmax=maximum(D)
    K = exp.(- D .^2 ./(Dmax^2*ϵ))
    
    k = Diagonal(vec(sum(K, dims=1)))
    L = inv(k)*K

    Λ=eigs(L, nev = maxoutdim, which = :LR)
    λ = real.(Λ[1])
    e = real.(Λ[2])
    Y = (λ) .* e'
    
    return L, λ, e, Y
end

function feature_vectors(states)
    Xs = [[p.x for p in S] for S in states]
    Vs = [[p.v for p in S] for S in states]
    # vmax = findmax([findmax(norm.(V))[1] for V in Vs])[1]

    features = momentConfig.(Xs,Vs)
    return Xs, Vs, features
end

function feature_vectors2(states)
    Xs = [[p.x for p in S] for S in states]
    Vs = [[p.v for p in S] for S in states]
    # vmax = findmax([findmax(norm.(V))[1] for V in Vs])[1]
    Xmeans = mean.(Xs)
    Xstds = std.(Xs)
    Vmeans = mean.(Vs)
    Vstds = std.(Vs)
    Zs = [[(Xs[i][j] .- Xmeans[i]) ./ Xstds[i] for j in eachindex(Xs[i])] for i in eachindex(Xs)]
    Ws = [[(Vs[i][j] .- Vmeans[i]) ./ Vstds[i] for j in eachindex(Vs[i])] for i in eachindex(Vs)]

    features = momentConfig.(Zs,Ws)
    return Xs, Vs, features
end

function distance_matrix(features)
    n = length(features)
    D = zeros(n,n)
    for i in eachindex(features)
        for j in 1:i-1
            D[i,j] = norm(features[i]-features[j])
            D[j,i] = D[i,j]
        end
    end
    return D
end


function SGECriterion(D,ϵ)
    Dmax=maximum(D)
   
    K1 = exp.(- D .^2 ./(Dmax^2*ϵ)) 
    k1 = Diagonal(vec(sum(K1, dims=1)))
    L1 = inv(k1)*K1

    K2 = exp.(- D .^2 ./(Dmax^2*(2*ϵ)))
    k2 = Diagonal(vec(sum(K2, dims=1)))
    L2 = inv(k2)*K2

    return sqrt(tr((L1*L1 - L2)*transpose(L1*L1 - L2)))
end

function SGECriterion2(D,ϵ)
    Dmax=maximum(D)
   
    K1 = exp.(- D .^2 ./(Dmax^2*ϵ)) 
    k1 = Diagonal(vec(sum(K1, dims=1)))
    L1 = inv(k1)*K1

    K2 = exp.(- D .^2 ./(Dmax^2*(2*ϵ)))
    k2 = Diagonal(vec(sum(K2, dims=1)))
    L2 = inv(k2)*K2

    return sqrt(real(eigs((L1*L1-L2)*transpose(L1*L1-L2),nev=1,which = :LR)[1][1]))
end

function SGECriteria(D,ϵs)
    SGE = [SGECriterion(D,ϵ) for ϵ in ϵs]
    maxSGE = maximum(SGE)
    return SGE ./ maxSGE
end

function SGECriteria2(D,ϵs)
    SGE = [SGECriterion2(D,ϵ) for ϵ in ϵs]
    maxSGE = maximum(SGE)
    return SGE ./ maxSGE
end

function ConnectivityCriterion(D,ϵ)
    Dmax=maximum(D)
    return length(findall(x->x<ϵ,D ./ Dmax))/length(D)
end

function ConnectivityCriteria(D,ϵs)
    return [ConnectivityCriterion(D,ϵ) for ϵ in ϵs]
end

function momentConfigVicsek(X,U)
    #Data set has X ∈ [0,L]×[0,L] and max(vᵢ) = v₀ = 3    
    # Z=[vcat(X[i]/10.0,U[i]/3.0) for i in 1:length(X)]
    Z=[vcat(X[i],U[i]) for i in eachindex(X)]
    meanZ = mean(Z)
    covZ = cov(Z)
    return [meanZ[1],meanZ[2],meanZ[3],meanZ[4],covZ[1,1],covZ[1,2],covZ[1,3],covZ[1,4],covZ[2,2],covZ[2,3],covZ[2,4],covZ[3,3],covZ[3,4],covZ[4,4]]
end