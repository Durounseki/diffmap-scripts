function designMatrix(X,F)
    return [f(x) for x in X, (key,f) in F]
end

function LSregresor(Y,Θ)
    β=(inv(transpose(Θ)*Θ)*transpose(Θ))*Y
    return β, (Y-Θ*β)⋅(Y-Θ*β)
end

multisets(m, k) = map(A -> [sum(A .== i) for i in 1:k], with_replacement_combinations(1:k, m)) #calculates the partitions of m in k integers allowing for 0 values

function CanonicalPolynomialBasis(n,m)
    # Return a dictionary of vector monomials of degrees 1,...,m
    F=OrderedDict{String,Function}("1"=>x->one(eltype(x)))
    for i in 1:m
        K=multisets(i,n)
        for k in K
            key=[@sprintf("x_1^%d ",k[1])]
            for j in 2:n
                push!(key,@sprintf("x_%d^%d ",j,k[j]))
            end
            F[join(key)]=X->.*((X.^k)...)
        end
    end
    return F
end

function TrigonometricBasis()
    F=OrderedDict{String,Function}()
    F["cos2"]=x->cos(x)^2
    F["sin2"]=x->sin(x)^2
    # F["1"]=x->one(eltype(x))
    # F["cos"]=x->cos(x)
    # F["sin"]=x->-sin(x)
    # F["cos2"]=x->0.5(3*cos(x)^2-1)
    # F["cossin"]=x->-3*cos(x)*sin(x)
    # F["sin2"]=x->3*sin(x)^2
    # F["cos3"]=x->0.5*(5*cos(x)^3-3*cos(x))
    # F["cos2sin"]=x->-1.5*(5*cos(x)^2-1)*sin(x)
    # F["cossin2"]=x->15*cos(x)*sin(x)^2
    # F["sin3"]=x->-15sin(x)^3
    return F
end

# function LegendreBasis()
#     F=OrderedDict{String,Function}()
