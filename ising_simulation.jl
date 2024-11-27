rand_ising2d(m, n=m) = rand(Int8[-1, 1], m, n)

Random.seed!(4649)

const β_crit = log(1+sqrt(2))/2
prob = @SVector [exp(-2*β_crit*k) for k in -4:4]

function ising2d_evolve!(s, β, niters; rng=default_rng())
    m, n = size(s)
    prob = @SVector [exp(-2*β*k) for k in -4:4]
    @fastmath @inbounds @simd for iter in 1:niters
        for j in 1:n 
            for i in 1:m
                let NN = s[ifelse(i == 1, m, i-1), j],
                    SS = s[ifelse(i == m, 1, i+1), j],
                    WW = s[i, ifelse(j == 1, n, j-1)],
                    EE = s[i, ifelse(j == n, 1, j+1)],
                    CT = s[i, j]
                    k = CT * (NN + SS + WW + EE)
                    s[i,j] = ifelse(rand(rng) < prob[k+5], -CT, CT)
                end
            end
        end
    end
end

function isingE(A,N;J=1)
    E=0
    for i in 0:N-1
        for j in 0:N-1
            E += -J*A[i+1,j+1]*(A[mod(i+1,N)+1,j+1]+A[mod(i-1,N)+1,j+1]+A[i+1,mod(j+1,N)+1]+A[i+1,mod(j-1,N)+1])#Nearest neighbours energy
        end
    end
    return E/2
end

function isingInitialEnsemble(M,N)
    return [rand_ising2d(N) for i in 1:M]
end

function isingDataSet(M,N,β,n_frames;n_t=1,rng=Random.default_rng())

    X=Vector{Matrix{Int64}}(undef,0)
    X0=isingInitialEnsemble(M,N)

    for i in 1:M
        push!(X,deepcopy(X0[i]))
        for j in 2:n_frames
            push!(X,deepcopy(X[end]))
            ising2d_evolve!(X[end],β,n_t;rng)
        end
    end

    return X
end

function isingDistanceMatrix(A)
    M = length(A)
    D=zeros(M,M)
    @sync for i in 1:M
        @spawn for j in 1:i-1
            D[i,j]=sum((A[i]-A[j]).^2)
            D[j,i]=D[i,j]
        end
    end
    return D
end

function isingDiffMapDistance(A,M,N,ϵ)
    
    # p=Progress(M*M)
    @sync for i in 1:M
        @spawn for j in 1:M
            D[i,j]=exp(-isingDistance(A[i],A[j])/((2*N)^2*ϵ))
            # next!(p)
        end
    end
end

function isingSteadyDataSet(M,N,β,n_frames,n_t;rng=default_rng())
    # X=[ising0condition(N)*(-1)^i for i in 1:M]
    X=isingInitialEnsemble(M,N)
    for x in X
        ising2d_evolve!(x,β,n_t;rng)
    end
    
    return X
end

function IsingSGEMinima(M,N,n_frames,Ts; n_t = 1, step = 10) #Fixed N data set for variable 
    minima = Vector{Vector{Int}}(undef,0)
    for T in Ts
        β = 1. / T
        X = isingDataSet(M,N,β,n_frames, n_t = n_t) #Data set at temperature T
        D = isingDistanceMatrix(X)
        ϵs = (0.00001:0.00001:0.01)∪(0.0101:0.0001:0.05)∪(0.051:0.001:0.1)∪(0.11:0.01:1)
        SGE = SGECriteria(D[1:step:end,1:step:end],ϵs)
        indices, minima_T = findminima(SGE)
    end
end

