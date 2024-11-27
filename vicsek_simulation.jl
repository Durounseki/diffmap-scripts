@inline function neighbors(d::Int)
    n = CartesianIndices((fill(-1:1, d)...,))
    return n[1:fld(length(n), 2)]
end

@inline function CellGrid(d::Int, L, r₀)
    indices = CartesianIndices((fill(-1:Int(fld(L,r₀)), d)...,))
    grid = Dict{CartesianIndex{d}, Vector{Int}}()
    for ci in indices
        grid[ci]=[]
    end
    return grid
end

function CellList(P, L, r₀, cells)
    n=length(P)
    d=length(P[1])
    for i in 1:n
        push!(cells[CartesianIndex(Int(fld(P[i][1],r₀)),Int(fld(P[i][2],r₀)))],i)
    end
    periodicBC(cells,L,r₀)
end

function ghostCells(L,r₀)
    tbBoundary=[CartesianIndex(i,Int(fld(L,r₀))-1) for i in 0:Int(fld(L,r₀))-1]
    append!(tbBoundary,[CartesianIndex(i,0) for i in 0:Int(fld(L,r₀))-1])
    rlBoundary=[CartesianIndex(Int(fld(L,r₀))-1,i) for i in 0:Int(fld(L,r₀))-1]
    append!(rlBoundary,[CartesianIndex(0,i) for i in 0:Int(fld(L,r₀))-1])
    return tbBoundary, rlBoundary
end

function periodicBC(cells,L,r₀)
    maxInd=Int(fld(L,r₀))
    tbBoundary, rlBoundary = ghostCells(L,r₀)
    for boundaryCell in tbBoundary
        if haskey(cells,boundaryCell)
            if boundaryCell[2]==0
                cells[boundaryCell+CartesianIndex(0,maxInd)] = cells[boundaryCell]
            elseif boundaryCell[2]==maxInd-1
                cells[boundaryCell+CartesianIndex(0,-maxInd)] = cells[boundaryCell]
            end
        end
    end
    for boundaryCell in rlBoundary
        if haskey(cells,boundaryCell)
            if boundaryCell[1]==0
                cells[boundaryCell+CartesianIndex(maxInd,0)] = cells[boundaryCell]
            elseif boundaryCell[1]==maxInd-1
                cells[boundaryCell+CartesianIndex(-maxInd,0)] = cells[boundaryCell]
            end
        end
    end
end

@inline function cellNeighbors!(ps, is, p, r)
    for (k, i) in enumerate(is[1:(end-1)])
        for j in is[(k+1):end]
            if norm(p[i]-p[j])<= r₀
                # println(i,j,p[i],p[j],norm(p[i]-p[j]))
                ps[i][j]=1
                ps[j][i]=1
            else
                ps[i][j]=0
                ps[j][i]=0
            end
        end
        ps[i][i]=1
    end
end

@inline function cellNeighbors!(ps, is, js, p, r)
    for i in is
        for j in js
            if norm(p[i]-p[j])<= r₀
                # println(i,j,p[i],p[j],norm(p[i]-p[j]))
                ps[i][j]=1
                ps[j][i]=1
            else
                ps[i][j]=0
                ps[j][i]=0
            end
        end
        ps[i][i]=1
    end
end

function near_neighbors(cells, P, r₀, FRN)
    offsets = neighbors(length(P[1]))
    # Iterate over non-empty cells
    for (cell, is) in cells
        # Pairs of points within the cell
        cellNeighbors!(FRN, is, P, r₀)
        # Pairs of points with non-empty neighboring cells
        for offset in offsets
            neigh_cell = cell + offset
            if haskey(cells, neigh_cell)
                @inbounds js = cells[neigh_cell]
                cellNeighbors!(FRN, is, js, P, r₀)
            end
        end
    end
end

function emptyCells(cells)
    for (ci,cell) in cells
        while length(cell)>0
            pop!(cell)
        end
    end
end

function fly(X,Θ,V,dt,L,r₀,N,v₀,η,cells,fixedRN)
    CellList(X, L, r₀, cells)
    # println("ok")
    near_neighbors(cells,X,r₀,fixedRN)
    # println("ok")
    
    Vmean=[normalize(mean(V[Is])) for Is in [findall(x->x>0,ff) for ff in fixedRN]]
    # println("ok")
    # Θnew=[u[2] > 0 ? acos(u[1])+η*(rand()-0.5) : 2π-acos(u[1])+η*(rand()-0.5) for u in Vmean]
    Θnew=[ atan(u[2],u[1])+η*(rand()-0.5) for u in Vmean]
    # println("ok")
    Vnew=[v₀*[cos(tt),sin(tt)] for tt in Θnew]
    # println("ok")
    
    for i in 1:N
        V[i]=Vnew[i]
        Θ[i]=Θnew[i]
        X[i]=X[i]+dt*Vnew[i]
    end
    # println("ok")
    PosPeriodicBC(X,L)
    # println("ok")
    
end

function PosPeriodicBC(X,L)
    for x in X
        x[1]=mod(x[1],L)
        x[2]=mod(x[2],L)
    end
end

function birdTrajectory(X,V,dt,N,nT,r₀,v₀,η,L,d)
    cells=CellGrid(d,L,r₀)
    fixedRN=[zeros(Int,N) for i in 1:N]
    filenameX=@sprintf("./vicsek_cells/positionsN%dL%2.lfr%2.lfv%2.lfeta%2.lf.csv",N,L,r₀,v₀,η)
    filenameV=@sprintf("./vicsek_cells/velocitiesN%dL%2.lfr%2.lfv%2.lfeta%2.lf.csv",N,L,r₀,v₀,η)
    Xio=open(filenameX,"w")
    Vio=open(filenameV,"w")
    writedlm(Xio, X, ", ")
    writedlm(Vio, V, ", ")
    close(Xio)
    close(Vio)
    Xio=open(filenameX,"a")
    Vio=open(filenameV,"a")
    for k in 1:nT
        fly(X,Θ,V,dt,L,r₀,N,v₀,η,cells,fixedRN)
        writedlm(Xio, X, ", ")
        writedlm(Vio, V, ", ")
        emptyCells(cells)
    end
    close(Xio)
    close(Vio)
end

function DataSet2(N,M,dt,nT,r₀,v₀,η,L,d)
    
    X₀=[[rand(2)*L for i in 1:N] for j in 1:M]
    Θ₀=[[2*π*rand() for i in 1:N] for j in 1:M]
    V₀=[[v₀*[cos(t),sin(t)] for t in tt] for tt in Θ₀]
    
    cells=CellGrid(d,L,r₀)
    fixedRN=[zeros(Int,N) for i in 1:N]
    for i in 1:N
        fixedRN[i][i]=1
    end
    
    for k in 1:nT
        fly(X₀[1],Θ₀[1],V₀[1],dt,L,r₀,N,v₀,η,cells,fixedRN)
        emptyCells(cells)
    end
    for i1 in 1:N
        for i2 in 1:N
            fixedRN[i1][i2]=0
        end
    end
    for i1 in 1:N
        fixedRN[i1][i1]=1
    end
    
    # S=[momentConfig(renormalisedData(X₀[1]),renormalisedData(V₀[1]))]
    S=[momentConfig(X₀[1] ./ L,V₀[1] ./ v₀)]
    
    for i in 2:M
        for k in 1:nT
            fly(X₀[i],Θ₀[i],V₀[i],dt,L,r₀,N,v₀,η,cells,fixedRN)
            emptyCells(cells)
        end
        
        for i1 in 1:N
        for i2 in 1:N
            fixedRN[i1][i2]=0
        end
        end
        for i1 in 1:N
            fixedRN[i1][i1]=1
        end
        # append!(S,[momentConfig(renormalisedData(X₀[i]),renormalisedData(V₀[i]))])
        append!(S,[momentConfig(X₀[i] ./ L,V₀[i] ./ v₀)])
    end

    return X₀, V₀, S
    
end

function dynamicDataset(N,M,dt,nT,r₀,v₀,η,L,d)
    X=Vector{Vector{Vector{Float64}}}(undef,0)
    X₀=[rand(2)*L for i in 1:N]
    Θ=Vector{Vector{Float64}}(undef,0)
    Θ₀=[2*π*rand() for i in 1:N]
    V=Vector{Vector{Vector{Float64}}}(undef,0)
    V₀=[v₀*[cos(t),sin(t)] for t in Θ₀]
    
    cells=CellGrid(d,L,r₀)
    fixedRN=[zeros(Int,N) for i in 1:N]
    for i in 1:N
        fixedRN[i][i]=1
    end
    
    push!(X,deepcopy(X₀))
    push!(Θ,deepcopy(Θ₀))
    push!(V,deepcopy(V₀))
    for k in 1:nT
        fly(X[end],Θ[end],V[end],dt,L,r₀,N,v₀,η,cells,fixedRN)
        emptyCells(cells)
    end
    for i1 in 1:N
        for i2 in 1:N
            fixedRN[i1][i2]=0
        end
    end
    for i1 in 1:N
        fixedRN[i1][i1]=1
    end
        
    S=[momentConfig((X[end] .- [[L/2,L/2]]) ./ L,V[end] ./ v₀)]
    
    for i in 2:M
        push!(X,deepcopy(X[end]))
        push!(Θ,deepcopy(Θ[end]))
        push!(V,deepcopy(V[end]))
        for k in 1:nT
            fly(X[end],Θ[end],V[end],dt,L,r₀,N,v₀,η,cells,fixedRN)
            emptyCells(cells)
        end
        for i1 in 1:N
            for i2 in 1:N
                fixedRN[i1][i2]=0
            end
        end
        for i1 in 1:N
            fixedRN[i1][i1]=1
        end
        append!(S,[momentConfig((X[end] .- [[L/2,L/2]]) ./ L,V[end] ./ v₀)])
        # S=[momentConfig(X[end] ./ L,V[end] ./ v₀)]
    end

    return X, V, S
    
end

function dynamicEnsemble(k,N,M,dt,nT,r₀,v₀,η,L,d)
    ensemble = [dynamicDataset(N,M,dt,nT,r₀,v₀,η,L,d)]
    for i in 2:k
        append!(ensemble,[dynamicDataset(N,M,dt,nT,r₀,v₀,η,L,d)])
    end
    return ensemble
end

#Minima for vicsek DataSet2 N=200;M=1000;dt=0.01;nT=200;r₀=2.0;v₀=3.0;η=0.05*π;L=10;d=2
#indices: 1083, 1432, 1473
#values: 0.0183, 0.082, 0.33

#Critical noise with N=200, L=10 is proportional to √(N/L²) = √2 ≃ 1.14
#But 0.5*π, π still orders
#I noticed that the simulation adds noise values between [-η/2,η/2] which means that the critical noise value we are going after is twice as large
#Therefore it makes sense that π still gives a phase transition

#η = 0.5π and the same values for the rest of the parameters is enough to generate a single trajectory that explores the entire slow manifold
# 1226
# 1482

function etaSet(N,M,dt,nT,r₀,v₀,ηs,L,d,ϵs,D;maxoutdim=3)
    
    filename1=@sprintf("./vicsekEvalsN%dL%2.lfr%.2lfv%.2lf.csv",N,L,r₀,v₀)
    Io1=open(filename1,"w")
    filename2=@sprintf("./vicsekOrderN%dL%2.lfr%.2lfv%.2lf.csv",N,L,r₀,v₀)
    Io2=open(filename2,"w")
    filename3=@sprintf("./vicsekSGEN%dL%2.lfr%.2lfv%.2lf.csv",N,L,r₀,v₀)
    Io3=open(filename3,"w")
       
    X,V,S=DataSet2(N,M,dt,nT,r₀,v₀,ηs[1],L,d)
    for i in eachindex(S)
        for j in 1:i-1
            D[i,j] = norm(S[i]-S[j])
            D[j,i] = D[i,j]
        end
    end
    SGE, ϵOpt = ϵOptimum(D,ϵs)
    transition_matrix, e_vals, e_vecs, modes = DiffMap(D, ϵOpt, maxoutdim = 3)
    writedlm(Io1, e_vals, ", ")
    close(Io1)
    writedlm(Io2, mean([sqrt(modes[2,i]^2+modes[3,i]^2) for i in 1:M]), ", ")
    close(Io2)
    writedlm(Io3, SGE, ", ")
    writedlm(Io3, ϵOpt)
    close(Io3)
    Io1=open(filename1,"a")
    Io2=open(filename2,"a")
    Io3=open(filename3,"a")
        
    for j in eachindex(ηs[2:end])
        X,V,S=DataSet2(N,M,dt,nT,r₀,v₀,ηs[j],L,d)
        for i in eachindex(S)
            for j in 1:i-1
                D[i,j] = norm(S[i]-S[j])
                D[j,i] = D[i,j]
            end
        end
        SGE, ϵOpt = ϵOptimum(D,ϵs)
        transition_matrix, e_vals, e_vecs, modes = DiffMap(D, ϵOpt, maxoutdim = maxoutdim)
        writedlm(Io1, e_vals, ", ")
        writedlm(Io2, mean([sqrt(modes[2,i]^2+modes[3,i]^2) for i in 1:M]), ", ")
        writedlm(Io3, SGE, ", ")
        writedlm(Io3, ϵOpt)
    end
    close(Io1)
    close(Io2)
    close(Io3)
end

function ϵOptimum(D,ϵs;ϵ_0=0.25)
    SGE = SGECriteria(D,ϵs)
    #Using the Peaks Library
    minima = findminima(SGE)
    if length(minima.indices) > 0
        return SGE, ϵs[minima.indices[end]];
    else
        return SGE, ϵ_0;
    end
end

function etaSet2(N,M,dt,nT,r₀,v₀,ηs,L,d,ϵs,ϵOpt,D;maxoutdim=3)
    
    filename1=@sprintf("./vicsekEvalsN%dL%2.lfr%.2lfv%.2lf.csv",N,L,r₀,v₀)
    Io1=open(filename1,"w")
    filename2=@sprintf("./vicsekOrderN%dL%2.lfr%.2lfv%.2lf.csv",N,L,r₀,v₀)
    Io2=open(filename2,"w")
       
    X,V,S=DataSet2(N,M,dt,nT,r₀,v₀,ηs[1],L,d)
    for i in eachindex(S)
        for j in 1:i-1
            D[i,j] = norm(S[i]-S[j])
            D[j,i] = D[i,j]
        end
    end
    
    transition_matrix, e_vals, e_vecs, modes = DiffMap(D, ϵOpt, maxoutdim = 3)
    writedlm(Io1, e_vals, ", ")
    close(Io1)
    writedlm(Io2, mean([sqrt(modes[2,i]^2+modes[3,i]^2) for i in 1:M]), ", ")
    close(Io2)
    Io1=open(filename1,"a")
    Io2=open(filename2,"a")
        
    for j in eachindex(ηs[2:end])
        X,V,S=DataSet2(N,M,dt,nT,r₀,v₀,ηs[j],L,d)
        for i in eachindex(S)
            for j in 1:i-1
                D[i,j] = norm(S[i]-S[j])
                D[j,i] = D[i,j]
            end
        end
        transition_matrix, e_vals, e_vecs, modes = DiffMap(D, ϵOpt, maxoutdim = maxoutdim)
        writedlm(Io1, e_vals, ", ")
        writedlm(Io2, mean([sqrt(modes[2,i]^2+modes[3,i]^2) for i in 1:M]), ", ")
    end
    close(Io1)
    close(Io2)
end