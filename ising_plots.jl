function plot_ϵCriteria(SGE,ϵs,iopt1, iopt2; w=400, h=350, ms = 10, label = false);

    F=Figure(size=(w,h))
    g = F[1,1:2] = GridLayout()

    ax = Axis(g[1, 1], aspect=1, yscale=log10)
    ylims!(ax,1e-3,10)
    xlims!(ax,0,1)

    lines!(ax,ϵs,SGE, color= :black, linewidth = 6)
    scatter!(ax,ϵs[iopt1],SGE[iopt1],marker = :dtriangle, color = Makie.wong_colors()[1], markersize=ms)
    scatter!(ax,ϵs[iopt2],SGE[iopt2],marker = :diamond, color = Makie.wong_colors()[2], markersize=ms)

    ax.xticklabelsvisible=label
    ax.yticklabelsvisible=label

    return F
end

function plot_ϵCriteria(SGE,ϵs,iopt; w=400, h=350, ms = 10, label = false);

    F=Figure(size=(w,h))
    g = F[1,1:2] = GridLayout()

    ax = Axis(g[1, 1], aspect=1, yscale = log10)
    ylims!(ax,1e-3,10)
    xlims!(ax,0,1)

    lines!(ax,ϵs,SGE, color = :black, linewidth = 6)
    
    scatter!(ax,ϵs[iopt],SGE[iopt],marker = :dtriangle, color = Makie.wong_colors()[1], markersize=ms)

    ax.xticklabelsvisible=label
    ax.yticklabelsvisible=label

    return F
end


function plot_projection(X,mode_ϵ_low,mode_ϵ_high; w=400, h=350, ms = 10, label = false)
    F=Figure(size=(w,h))
    g = F[1,1:2] = GridLayout()

    ax = Axis(g[1, 1], aspect=1)

    scatter!(ax, mode_ϵ_low, mean.(X), color=Makie.wong_colors()[1], marker = :dtriangle, markersize=ms)
    scatter!(ax, mode_ϵ_high, mean.(X), color=Makie.wong_colors()[2], marker = :diamond, markersize=ms)

    theta=designMatrix(mode_ϵ_high,CanonicalPolynomialBasis(1,1))
    beta=LSregresor(mean.(X),theta)

    js = 1:length(mode_ϵ_high)

    lines!(ax,mode_ϵ_high, mode_ϵ_high .* beta[1][2] .+ beta[1][1], color = :black, linewidth = 6)

    ax.xticklabelsvisible=label
    ax.yticklabelsvisible=label

    return F

end

function plot_projection2(X,mode_ϵ_low,mode_ϵ_high; w=400, h=350, ms = 10, label = false)
    F=Figure(size=(w,h))
    g = F[1,1:2] = GridLayout()

    ax = Axis(g[1, 1], aspect=1)

    scatter!(ax, mode_ϵ_low, X, color=Makie.wong_colors()[1], marker = :dtriangle, markersize=ms)
    scatter!(ax, mode_ϵ_high, X, color=Makie.wong_colors()[2], marker = :diamond, markersize=ms)

    theta=designMatrix(mode_ϵ_low,CanonicalPolynomialBasis(1,1))
    beta=LSregresor(mean.(X),theta)

    lines!(ax,-0.01:0.001:0.01, x -> x * beta[1][2] + beta[1][1], color = :black, linewidth = 6)

    theta=designMatrix(mode_ϵ_high,CanonicalPolynomialBasis(1,1))
    beta=LSregresor(mean.(X),theta)

    lines!(ax,-0.01:0.001:0.01,x-> beta[1][2] * x + beta[1][1], linestyle = :dash, dashpattern =(10,30), color = :black, linewidth = 6)

    ax.xticklabelsvisible=label
    ax.yticklabelsvisible=label

    return F

end

function plot_projection(X,mode; w=400, h=350, ms = 10, label = false)
    F=Figure(size=(w,h))
    g = F[1,1:2] = GridLayout()

    ax = Axis(g[1, 1], aspect=1)

    scatter!(ax, mode, mean.(X), color=Makie.wong_colors()[1], marker = :dtriangle, markersize=ms)

    theta=designMatrix(mode,CanonicalPolynomialBasis(1,1))
    beta=LSregresor(mean.(X),theta)

    lines!(ax,mode, mode .* beta[1][2] .+ beta[1][1], color = :black, linewidth = 6)

    ax.xticklabelsvisible=label
    ax.yticklabelsvisible=label

    return F

end

function plot_projection(mode_ϵ_low_2,mode_ϵ_low_3,mode_ϵ_hi_2,mode_ϵ_hi_3,order_param; w=400, h=350, ms = 10, label = false)
    F=Figure(size=(w,h))
    g = F[1,1:2] = GridLayout()

    ax = Axis(g[1, 1], aspect=1)

    colormap1 = RGBAf.(Colors.color.(to_colormap(:roma)), 0.25)

    scatter!(ax, mode_ϵ_low_2, mode_ϵ_low_3, color=order_param, marker = :dtriangle, markersize=ms, colormap=colormap1)
    scatter!(ax, mode_ϵ_hi_2, mode_ϵ_hi_3, color=order_param, marker = :diamond, markersize=ms, colormap=colormap1)

    # theta=designMatrix(mode_ϵ_high,CanonicalPolynomialBasis(1,1))
    # beta=LSregresor(mean.(X),theta)

    # js = 1:length(mode_ϵ_high)

    # lines!(ax,mode_ϵ_high, mode_ϵ_high .* beta[1][2] .+ beta[1][1], color = :black, linewidth = 6)

    ax.xticklabelsvisible=label
    ax.yticklabelsvisible=label

    return F

end

function plot_projection3(mode_2,mode_3,order_param; w=400, h=350, ms = 10, label = false)
    F=Figure(size=(w,h))
    g = F[1,1:2] = GridLayout()

    ax = Axis(g[1, 1], aspect=1)

    colormap1 = RGBAf.(Colors.color.(to_colormap(:roma)), 0.75)

    scatter!(ax, mode_2, mode_3, color=order_param, marker = :dtriangle, markersize=ms, colormap=colormap1)
    # scatter!(ax, mode_ϵ_hi_2, mode_ϵ_hi_3, color=order_param, marker = :diamond, markersize=ms, colormap=colormap1)

    # theta=designMatrix(mode_ϵ_high,CanonicalPolynomialBasis(1,1))
    # beta=LSregresor(mean.(X),theta)

    # js = 1:length(mode_ϵ_high)

    # lines!(ax,mode_ϵ_high, mode_ϵ_high .* beta[1][2] .+ beta[1][1], color = :black, linewidth = 6)

    ax.xticklabelsvisible=label
    ax.yticklabelsvisible=label

    return F

end

function plot_series(X,mode_ϵ_low,mode_ϵ_high; w=400, h=350, ms = 10, label = false)
    F=Figure(size=(w,h))
    g = F[1,1:2] = GridLayout()

    ax = Axis(g[1, 1], aspect=1)
    js = 1:length(X)

    lines!(ax, js, mean.(X), color=:black, linewidth = 4)
    
    theta=designMatrix(mode_ϵ_high,CanonicalPolynomialBasis(1,1))
    beta=LSregresor(mean.(X),theta)

    scatter!(ax, js[1:20:end], mode_ϵ_high[1:20:end] .* beta[1][2] .+ beta[1][1], color=Makie.wong_colors()[2], marker = :diamond, markersize=ms)
    scatter!(ax, js[1:20:end], mode_ϵ_low[1:20:end] * sqrt(length(js)), color=Makie.wong_colors()[1], marker = :dtriangle, markersize=ms)

    ax.xticklabelsvisible=label
    ax.yticklabelsvisible=label

    return F

end

function plot_series2(X,mode_ϵ_low,mode_ϵ_high; w=400, h=350, ms = 10, label = false)
    F=Figure(size=(w,h))
    g = F[1,1:2] = GridLayout()

    ax = Axis(g[1, 1], aspect=1)
    js = 1:length(X)

    lines!(ax, js, mean.(X), color=:black, linewidth = 4)
    
    theta=designMatrix(mode_ϵ_low,CanonicalPolynomialBasis(1,1))
    beta=LSregresor(X,theta)
    scatter!(ax, js[1:2:end], mode_ϵ_low[1:2:end] .* beta[1][2] .+ beta[1][1], color=Makie.wong_colors()[1], marker = :dtriangle, markersize=ms)

    theta=designMatrix(mode_ϵ_high,CanonicalPolynomialBasis(1,1))
    beta=LSregresor(X,theta)
    scatter!(ax, js[1:2:end], mode_ϵ_high[1:2:end] .* beta[1][2] .+ beta[1][1], color=Makie.wong_colors()[2], marker = :diamond, markersize=ms)

    ax.xticklabelsvisible=label
    ax.yticklabelsvisible=label

    return F
end

function plot_series2(X,mode; w=400, h=350, ms = 10, label = false)
    F=Figure(size=(w,h))
    g = F[1,1:2] = GridLayout()

    ax = Axis(g[1, 1], aspect=1)
    js = 1:length(X)
    ylims!(ax,-1.05,1.05)

    lines!(ax, js, mean.(X), color=:black, linewidth = 2)
    
    theta=designMatrix(mode,CanonicalPolynomialBasis(1,1))
    beta=LSregresor(X,theta)
    scatter!(ax, js[1:2:end], mode[1:2:end] .* beta[1][2] .+ beta[1][1], color=Makie.wong_colors()[1], marker = :dtriangle, markersize=ms)

    ax.xticklabelsvisible=label
    ax.yticklabelsvisible=label

    return F
end

function plot_seriesNorm(X,mode_ϵ_low,mode_ϵ_high; w=400, h=350, ms = 10, label = false)
    F=Figure(size=(w,h))
    g = F[1,1:2] = GridLayout()

    ax = Axis(g[1, 1], aspect=1)
    js = 1:length(X)

    lines!(ax, js, mean.(X), color=:black, linewidth = 4)
    
    x1min=findmin(mode_ϵ_low)[1]; x1max=findmax(mode_ϵ_low)[1]
    x1 = (2*(mode_ϵ_low) .- (x1max+x1min)) ./ (x1max-x1min);
    x2min=findmin(mode_ϵ_high)[1]; x2max=findmax(mode_ϵ_high)[1]
    x2 = (2*(mode_ϵ_high) .- (x2max+x2min)) ./ (x2max-x2min);

    scatter!(ax, js[1:10:end], x1[1:10:end], color=Makie.wong_colors()[1], marker = :dtriangle, markersize=ms)
    scatter!(ax, js[1:10:end], x2[1:10:end], color=Makie.wong_colors()[2], marker = :diamond, markersize=ms)

    ax.xticklabelsvisible=label
    ax.yticklabelsvisible=label

    return F
end

function plot_series(X,mode; w=400, h=350, ms = 10, label = false)
    F=Figure(size=(w,h))
    g = F[1,1:2] = GridLayout()

    ax = Axis(g[1, 1], aspect=1)
    js = 1:length(X)
    ylims!(ax,-1.05,1.05)
    
    theta=designMatrix(mode,CanonicalPolynomialBasis(1,1))
    beta=LSregresor(mean.(X),theta)

    lines!(ax, js, mean.(X), color=:black, linewidth = 4)
    scatter!(ax, js, mode .* beta[1][2] .+ beta[1][1], color=Makie.wong_colors()[1], marker = :dtriangle, markersize=ms)
    
    # lines!(ax, js, mode .* beta[1][2] .+ beta[1][1], color=Makie.wong_colors()[1], marker = :dtriangle, markersize=ms)

    ax.xticklabelsvisible=label
    ax.yticklabelsvisible=label

    return F

end

function plot_spectrum(evals_ϵ_low,evals_ϵ_high;w=400, h=350, ms = 10, label = false)
    F=Figure(size=(w,h))
    g = F[1,1:2] = GridLayout()

    ax = Axis(g[1, 1], aspect=1)
    scatter!(ax,evals_ϵ_low[1:10],marker = :dtriangle, color = Makie.wong_colors()[1], markersize=ms)
    scatter!(ax,evals_ϵ_high[1:10],marker = :diamond, color = Makie.wong_colors()[2], markersize=ms)

    ax.xticklabelsvisible=label
    ax.yticklabelsvisible=label

    return F
end

function plot_spectrum(evals;w=400, h=350, ms = 10, label = false)
    F=Figure(size=(w,h))
    g = F[1,1:2] = GridLayout()

    ax = Axis(g[1, 1], aspect=1)
    scatter!(ax,evals[1:10],marker = :dtriangle, color = Makie.wong_colors()[1], markersize=ms)

    ax.xticklabelsvisible=label
    ax.yticklabelsvisible=label

    return F
end