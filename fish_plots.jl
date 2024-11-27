function BSplineAppex(X)
    n=length(X)
    
    c1=[1.0/4]
    S=[(6*X[2]-X[1])/4]
    for i in 2:n-3
        append!(S,[(6*X[i+1]-S[end])/(4-c1[end])])
        append!(c1,1/(4-c1[end]))
    end
    append!(S,[(6*X[n-1]-X[n]-S[end])/(4-c1[end])])
    append!(c1,0)
    
    B=[S[end]]
    for i in 1:length(S)-1
        append!(B,[S[end-i]-c1[end-i]*B[end]])
    end 
    return [B[end-i] for i in 0:length(B)-1]
end

function lenIs(I)
    if length(I)>0
        D=[[I[1]-1,1]]
        for i in I[2:end]
            if i-(sum(D[end]))==1
                D[end][2]+=1
            else
                append!(D,[[i-1,1]])
            end
        end
        return D
    else 
        return [[0,0]]
    end
end

function BSCubicInterpolator(X)
    B=BSplineAppex(X)
    Q=[X[1]]
    n=length(X)
    append!(Q,[X[1]+(B[1]-X[1])/3])
    append!(Q,[X[1]+2*(B[1]-X[1])/3])
    append!(Q,[X[2]])
    append!(Q,[X[2]])
    for i in 2:n-2
        append!(Q,[B[i-1]+(B[i]-B[i-1])/3])
        append!(Q,[B[i-1]+2*(B[i]-B[i-1])/3])
        append!(Q,[X[i+1]])
        append!(Q,[X[i+1]])
    end
    append!(Q,[B[n-2]+(X[n]-B[n-2])/3])
    append!(Q,[B[n-2]+2*(X[n]-B[n-2])/3])
    append!(Q,[X[n]])
    return Q
end

function trajectoryInterpolator(X)
    Is1=findall(x->norm(x)==0,X)
    BS=BSCubicInterpolator(X[setdiff(1:1:length(X),Is1)])
    V=[3*(BS[(i-1)*4+2]-BS[(i-1)*4+1]) for i in 1:length(setdiff(1:1:length(X),Is1))-1]
    append!(V,[3*(BS[end]-BS[end-1])])
    Is2=lenIs(Is1)
    s=0
    for i in eachindex(Is2)
        
        if Is2[i][1]==0
            for j in 1:Is2[i][2]
                X[Is2[i][1]+(Is2[i][2]-j+1)]=2*X[Is2[i][1]+(Is2[i][2]-j+1)+1]-X[Is2[i][1]+(Is2[i][2]-j+1)+2]
                insert!(V,1,V[1])
            end
        elseif sum(Is2[i])==length(X)
            for j in 1:Is2[i][2]
                X[Is2[i][1]+j]=2*X[Is2[i][1]+j-1]-X[Is2[i][1]+j-2]
                insert!(V,length(V)+1,V[end])
            end
        else
            h=1.0/(Is2[i][2]+1)
            k=Is2[i][1]-s
            for j in 1:Is2[i][2]
                X[Is2[i][1]+j]=(1-j*h)^3*BS[(k-1)*4+1]+3*(1-j*h)^2*j*h*BS[(k-1)*4+2]+3*(1-j*h)*(j*h)^2*BS[(k-1)*4+3]+(j*h)^3*BS[(k-1)*4+4]
                insert!(V,Is2[i][1]+j,3*(1-j*h)^2*(BS[(k-1)*4+2]-BS[(k-1)*4+1])+6*(1-j*h)*j*h*(BS[(k-1)*4+3]-BS[(k-1)*4+2])+3*(j*h)^2*(BS[(k-1)*4+4]-BS[(k-1)*4+3]))
            end
        end
        s+=Is2[i][2]
    end
    return V
end

function fishPositions(N,m)
    filenameX=@sprintf("./fish/trajectories%d_%02d.csv",N,m)
    Xio=open(filenameX,"r")
    
    X=[]
    for line in eachline(Xio)
        append!(X,[parse.(Float64,split(line,","))])
    end
    close(Xio)
    
    Y=[[X[(i-1)*N+j] for i in 1:Int(length(X)/N)] for j in 1:N]
    
    V=[trajectoryInterpolator(y) for y in Y]
    return [[Y[j][i] for j in 1:N] for i in 1:length(Y[1])], [[V[j][i] for j in 1:N] for i in 1:length(V[1])]
end

function renormalisedData(X)
    xMin=findmin([x[1] for x in X])[1]
    yMin=findmin([x[2] for x in X])[1]
    xMax=findmax([x[1] for x in X])[1]
    yMax=findmax([x[2] for x in X])[1]
    return [[(2*x[1]-(xMin+xMax))/(xMax-xMin),(2*x[2]-(yMax+yMin))/(yMax-yMin)] for x in X]
end

function fishFeatures(N,m)

    X=fishPositions(N,m)[1]

    V=noisyDerivative(X)
    X_2=X[4:end-3]

    Xn=[renormalisedData(x) for x in X_2]
    Vn=[renormalisedData(x) for x in V]

    MXV=[momentConfig(Xn[i],Vn[i]) for i in eachindex(Xn)]

    return X_2, V, Xn, Vn, MXV
    
end

function fishAccelerationFeatures(N,m)

    X=fishPositions(N,m)[1]

    V=noisyDerivative(X)
    A=noisyDerivative(V)
    X_2=X[7:end-6]
    V_2=V[4:end-3]

    Xn=[renormalisedData(x) for x in X_2]
    Vn=[renormalisedData(v) for v in V_2]
    An=[renormalisedData(a) for a in A]

    MXVA = [AccMomentConfig(Xn[i],Vn[i],An[i]) for i in eachindex(Xn)]

    return X_2, V_2, A, Xn, Vn, An, MXVA

end

function fishZscoreFeatures(N,m; acc=false)
    
    X=fishPositions(N,m)[1]

    V=noisyDerivative(X)
    X_1=X[4:end-3]
    
    if acc

        A=noisyDerivative(V)
        X_2=X[7:end-6]
        V_2=V[4:end-3]
        MXVA = [zScoreFeatures(X_2[i],V_2[i],A[i]) for i in eachindex(A)]
        return X_2, V_2, A, MXVA

    else

        MXV = [zScoreFeatures(X_1[i],V[i]) for i in eachindex(V)]
        return X_1, V, MXV

    end

end

function fishNN(X)
    return [findmin([norm(x1-x2) for x2 in X]) for x1 in X]
end

function fishRelOrientation(V)
    return [v1⋅mean([v2 for v2 in V]) for v1 in V]
end

function spectralGaps(e_vals)
    return [e_vals[i]-e_vals[i+1] for i in eachindex(e_vals[1:end-1])]
end

function spectralGapCriteria(D,ϵs;maxoutdim=10)
    SGϵ = [[] for i in 1:maxoutdim]
    for ϵ in ϵs
        transition_matrix, e_vals, e_vecs, modes = DiffMap(D,ϵ,maxoutdim)
        SGs = spectralGaps(e_vals)
        for i in eachindex(SGs)
            push!(SGϵ[i],SGs[i])
        end
    end
    return SGϵ
end

function noisyDerivative(X)
    return [(5*(X[i+1]-X[i-1])+4*(X[i+2]-X[i-2])+X[i+3]-X[i-3])/(2^5) for i in 4:length(X)-3]
end

function triangle(x,v,s)
    Vp=normalize([-v[2],v[1]])
    Vn=normalize(v)
    return Polygon(Point2f[(x[1]-0.5*s*Vp[1],x[2]-0.5*s*Vp[2]),(x[1]+0.5*s*Vp[1],x[2]+0.5*s*Vp[2]),(x[1]+2*s*Vn[1],x[2]+2*s*Vn[2])])
end

function fishSnapshot(X,V,time;scale=0.05, shadow=10, width=400, height=400,kwargs...)
    
    school_config = triangle.(X[time],V[time],scale)
    school_shadows = [triangle.(X[time-i],V[time-i],scale) for i in 1:shadow-1]
    shadows_color = [RGBAf(1-(shadow-i)/shadow, 1-(shadow-i)/shadow, 1-(shadow-i)/shadow,(shadow-i)/shadow) for i in 0:shadow-1]
    
    F = Figure(size=(width,height), aspectratio=1)
    ax=Axis(F[1,1])
    for i in 2:shadow
    #     # school_config = triangle.(X[time-shadow],V[time-shadow],scale)
        poly!(school_shadows[shadow+1-i],color=shadows_color[shadow+1-i])
        # poly!(school_config,color=shadows_color[2])
    end
    
    poly!(ax,school_config,color=:black)
    lines!([sqrt(2) .* (cos(q),sin(q)) for q in 0:0.1:(2*π+0.1)],color=:black,linewidth=3)
    hidedecorations!(ax)

    return F
    # return shadows_color
end

function fishSnapshots(N,m,times;kwargs...)
    X, V, Xn, Vn, features = fishFeatures(N,m)
    for t in times
        filename=@sprintf("./fish/snapshotN%dm%dt%d.pdf",N,m,t)
        F=fishSnapshot(Xn,Vn,t;kwargs...)
        save(filename,F)
    end
end

function fishDiffmapPlot(e_vals,ϵs,ϵopt,Rs,Rmin;w=800, h=350, ms = 30, legend = true,axDec = true)
    F=Figure(size=(800,800))
    g1=F[1,1:2]=GridLayout()
    g2=F[1,3:4]=GridLayout()
    axg1 = Axis(g1[1,1],aspect=1)
    axg2 = Axis(g2[1,1],aspect=1)

    Tags=[L"N=100",L"N=80",L"N=60"]
    markers=[:circle,:rect,:utriangle]
    linetypes=[nothing,:dash,:dot]
    
    for i in eachindex(Rs)
        scatter!(axg1,1:10,e_vals[i][1:10],color=Makie.wong_colors()[i], marker=markers[i], markersize=ms, legend=false)
        if legend
            scatter!(axg1,12,0,marker=markers[i],color=color=Makie.wong_colors()[i], markersize = ms/2, label=Tags[i])
        end
        lines!(axg2,ϵs,Rs[i],linewidth = 6, linestyle=linetypes[i])
    end
    for i in eachindex(Rs)
        scatter!(axg2,ϵopt[i],Rmin[i],marker=markers[i],markersize=2*ms/3)
    end
    if legend
        xlims!(axg1,0,10.2)
        axislegend(axg1)
    end
    axg1.xticklabelsvisible=axDec
    axg1.yticklabelsvisible=axDec
    axg2.xticklabelsvisible=axDec
    axg2.yticklabelsvisible=axDec
    return F

end

function fishTimeseriesPlot(features,order;ms=5,axDec = true)
    F=Figure(size=(800,800))
    g1=F[1,1:2]=GridLayout()
    g2=F[1,3:4]=GridLayout()
    axg1 = Axis(g1[1,1],aspect=1)
    axg2 = Axis(g2[1,1],aspect=1)



    L = [f[8]-f[10] for f in features]
    theta=designMatrix(order,CanonicalPolynomialBasis(1,1))
    beta=LSregresor(L,theta)
    ts = (1:length(L))/60

    scatter!(axg1,order,L, markersize=ms)
    lines!(axg1,order, order .* beta[1][2] .+ beta[1][1], color = :black, linewidth = 6)

    lines!(axg2,ts,L,linewidth=6)
    lines!(axg2,ts[1:10:end],order[1:10:end] .* beta[1][2] .+ beta[1][1],linestyle=:dash,linewidth=3)
    
    
    axg1.xticklabelsvisible=axDec
    axg1.yticklabelsvisible=axDec

    axg2.xticklabelsvisible=axDec
    axg2.yticklabelsvisible=axDec

    return F

end


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

#Fish diffmap02 has n=60 m=1
#Fish diffmap03 has n=60 m=2
#Fish diffmap04 has n=60 m=3