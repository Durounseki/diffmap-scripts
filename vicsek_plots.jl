# function plot_ϵCriteria(SGE,ϵs,iopt1, iopt2; w=400, h=350, ms = 10, label = false);

#     F=Figure(size=(w,h))
#     g = F[1,1:2] = GridLayout()

#     ax = Axis(g[1, 1], aspect=1, yscale=log10)
#     ylims!(ax,1e-5,10)
#     xlims!(ax,0,1)

#     lines!(ax,ϵs,SGE,linewidth = 2)
#     scatter!(ax,ϵs[iopt1],SGE[iopt1],marker = :dtriangle, color = Makie.wong_colors()[1], markersize=ms)
#     scatter!(ax,ϵs[iopt2],SGE[iopt2],marker = :diamond, color = Makie.wong_colors()[2], markersize=ms)

#     ax.xticklabelsvisible=label
#     ax.yticklabelsvisible=label

#     return F
# end

# function plot_ϵCriteria(SGE,ϵs,iopt; w=400, h=350, ms = 10, label = false);

#     F=Figure(size=(w,h))
#     g = F[1,1:2] = GridLayout()

#     ax = Axis(g[1, 1], aspect=1, yscale = log10)
#     ylims!(ax,1e-5,10)
#     xlims!(ax,0,1)

#     lines!(ax,ϵs,SGE, color = :black, linewidth = 6)
    
#     scatter!(ax,ϵs[iopt],SGE[iopt],marker = :dtriangle, color = Makie.wong_colors()[1], markersize=ms)

#     ax.xticklabelsvisible=label
#     ax.yticklabelsvisible=label

#     return F
# end

# function plot_order(vx,vy,mode2,mode3; w=400, h=350, ms = 10, label=false)
#     F=Figure(size=(800,800))
#     g1=F[1,1:2]=GridLayout()
#     g2=F[1,3:4]=GridLayout()
#     axg1 = Axis(g1[1,1],aspect=1)
#     axg2 = Axis(g2[1,1],aspect=1)

#     cosv = [([vx[i],vy[i]]⋅[vx[1],vy[1]])/(norm([vx[i],vy[i]])*norm([vx[1],vy[1]])) for i in eachindex(vx)]
#     cosψ = [[mode2[i],mode3[i]]⋅[mode2[1],mode3[1]]/(norm([mode2[i],mode3[i]])*norm([mode2[1],mode3[1]])) for i in eachindex(mode2)]

#     theta=designMatrix(cosψ,CanonicalPolynomialBasis(1,1))
#     beta=LSregresor(cosv,theta)

#     js = 1:length(cosv)

#     scatter!(axg1,cosψ,cosv, markersize=ms)
#     # lines!(axg1,cosψ, cosψ .* beta[1][2] .+ beta[1][1], color = :black, linewidth = 6)

#     axg1.xticklabelsvisible=label
#     axg1.yticklabelsvisible=label

#     lines!(axg2,js,cosv,linewidth=6)
#     lines!(axg2,js,cosψ,linestyle=:dash,linewidth=3)
#     # lines!(axg2,js[1:10:end],cosψ[1:10:end] .* beta[1][2] .+ beta[1][1],linestyle=:dash,linewidth=3)

#     axg2.xticklabelsvisible=label
#     axg2.yticklabelsvisible=label

#     return F
# end

# function generate_features(Xfilename,Vfilename)
#     #The data consists of the positions and velocities of N=200 particles over M=1000 snapshots, taken every 2s
#     positions_matrix = readdlm(Xfilename,',')
#     velocities_matrix = readdlm(Vfilename,',')
#     N = 200; M = 1000
#     #The data is stored as a matrix of size (NM) x 2
#     positions_array = [[[positions_matrix[i+(j-1)*N,1],positions_matrix[i+(j-1)*N,2]] for i in 1:N] for j in 1:M]
#     velocities_array = [[[velocities_matrix[i+(j-1)*N,1],velocities_matrix[i+(j-1)*N,2]] for i in 1:N] for j in 1:M]
#     features = [momentConfigVicsek(positions_array[i],velocities_array[i]) for i in 1:M]
#     return features
# end
    
# function read_features(filename)
#     features = readdlm(filename,',')
#     for i in eachindex(features[:,1])
#         features[i,1] = parse(Float64,features[i,1][2:end])
#         features[i,14] = parse(Float64,features[i,14][1:end-1])
#         features[i,15] = parse(Float64,features[i,15][2:end])
#         features[i,end] = parse(Float64,features[i,end][1:end-1])
#     end
#     f1 = [[f for f in features[i,:][1:14]] for i in eachindex(features[:,1])]
#     f2 = [[f for f in features[i,:][15:end]] for i in eachindex(features[:,1])]
#     f = vcat(f1,f2)
#     return f
# end

# #N=200;M=10000;dt=0.01;nT=200;r₀=2.0;v₀=3.0;η=π/6;L=10;d=2
# #X4, V4, S4 = dynamicDataset(N,M,dt,nT,r₀,v₀,η,L,d)
# #Ellipse fitting
# function fitzgibbon_ellipse(x,y)
#     D = [x.^2 x.*y y.^2 x y ones(length(x))]
#     S = D' * D
#     C = zeros(6,6)
#     C[1,3] = C[3,1] = 2
#     C[2,2] = -1
#     eigvals, eigvecs = eigen(S,C)
#     coeffs = eigvecs[:,argmin(eigvals)]
#     coeffs /= sqrt(4*coeffs[1]*coeffs[3]-coeffs[2]^2);
#     return D, coeffs
# end

function plot_transition_fixedL(Ns,L,ηs,r₀,v₀;w=800,h=350, ms=10, label=false)
    markers=[:circle, :rect, :diamond, :hexagon, :utriangle, :cross, :star, :pentagon]
    F = Figure(size = (w,h))
    g = F[1,1:2] = GridLayout()

    ax = Axis(g[1,1], aspect=2)
    # xlims!(ax,(0,7))

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
    # ax2 = Axis(g[1, 1], xaxisposition=:top, xticks=(jumps,[@sprintf("%.2lf",jump) for jump in jumps]), aspect=2)
    # linkxaxes!(ax,ax2)
    # hideydecorations!(ax2, ticks = false)
    # hideydecorations!(ax2, grid = false)
    
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