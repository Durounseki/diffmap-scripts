function plot_transition_fixedL(Ns,L,ηs,r₀,v₀;w=800,h=350, ms=10, label=false)
    markers=[:circle, :rect, :diamond, :hexagon, :utriangle, :cross, :star, :pentagon]
    F = Figure(size = (w,h))
    g = F[1,1:2] = GridLayout()

    ax = Axis(g[1,1], aspect=2)

    legend_entries = []
    jumps = []
    for i in eachindex(Ns)
        filename = @sprintf("./vicsekOrderN%dL%2.lfr%.2lfv%.2lf.csv",Ns[i],L,r₀,v₀)
        Io = open(filename,"r")
        order=readdlm(Io,',',Float64,'\n')
        treshold = (findmax(order)[1]-findmin(order)[1])/4
        jump_index = findSingleJump(order,treshold)
        push!(jumps,ηs[jump_index])
        scatter_obj = scatter!(ax,ηs,order[1:end],marker=markers[(i-1)%length(markers) + 1],markersize=ms)
        vlines!(ax,jumps[end],linestyle=:dash,linewidth=3, color=:black)
        push!(legend_entries,scatter_obj)
    end
    
    labels=[@sprintf("N=%d",N) for N in Ns]
    Legend(F[1,3],legend_entries,labels)

    ax.xticklabelsvisible=label
    ax.yticklabelsvisible=label

    return F
end

function plot_scaling_fixedL(Ns,L,ηs,r₀,v₀;w=800,h=350, ms=10, label=false)
    F = Figure(size = (w,h))
    g = F[1,1:2] = GridLayout()

    ax = Axis(g[1,1], aspect=1, yscale=log10, xscale=log10,yminorticksvisible = true,xminorticksvisible = true, yminorticks = IntervalsBetween(3),xminorticks = IntervalsBetween(10))
    xlims!(ax,(50,2000))
    ylims!(ax,(π/4,2π+1))
    ax.xticks=[100,400,1600]
    ax.yticks=([π/4,π/2,π,2*π],["π/4","π/2","π","2π"])
    jumps = []
    for i in eachindex(Ns)
        filename = @sprintf("./vicsekOrderN%dL%2.lfr%.2lfv%.2lf.csv",Ns[i],L,r₀,v₀)
        Io = open(filename,"r")
        order=readdlm(Io,',',Float64,'\n')
        treshold = (findmax(order)[1]-findmin(order)[1])/4
        jump_index = findSingleJump(order,treshold)
        push!(jumps,ηs[jump_index])
    end
    
    scatter!(ax,Ns,[jump for jump in jumps],markersize=ms)
    
    return F

end

function findSingleJump(X,treshold)
    for i in eachindex(X[1:end-1])
        if abs(X[i+1] - X[i]) > treshold
            return i
        end
    end
end