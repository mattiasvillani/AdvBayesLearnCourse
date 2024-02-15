using Plots
using Distributions
using LaTeXStrings
using ColorSchemes
import ColorSchemes: Paired_12; colors = Paired_12
basicColors = Paired_12[[1,2,7,8,3,4,5,6,9,10]]

default(legend = true)
figFolder = "/home/mv/Dropbox/Teaching/AdvBayesLearnCourse/Slides/Figures/"

function logmarglike(xBar, n, σ, κ0)
        return -(1/2)*log(κ0/(κ0 + n)) - ((n*xBar^2)/(2*σ^2))*(n/(κ0 + n))
end

# Data
xBar = 0.6
n = 10
σ = 1
μ0 = 0

κ0s = 0.01:0.01:1000
BFs = [exp(logmarglike(xBar, n, σ, κ0)) for κ0 in κ0s]
plot(κs, BFs, lw = 3, ylabel = "Bayes factor", xlabel = L"\kappa_0", 
        ylims = [0,5], xaxis=:log, label = "")
savefig(figFolder*"NormalMargLikeAllKappa0.pdf")


θGrid = -3:0.01:3
margDensH0 = pdf.(Normal(0,σ*√(1/n)), θGrid)
p = plot(xlabel = L"\theta", ylabel = "Marginal likelihood")
plot!(p, θGrid, margDensH0, color = basicColors[6], lw = 2, label = L"H_0")
for (count, κ0) ∈ enumerate([100 10 1 0.1])
        margDensH1 = pdf.(Normal(μ0,σ*√(1/n + 1/κ0)), θGrid)
        plot!(p, θGrid, margDensH1, color = basicColors[count], lw = 2,
                label = LaTeXString("\$\\kappa_0 = $κ0\$"))
end
scatter!([xBar], [0], seriestype = :scatter, color = basicColors[10], 
        label = L"\bar x", markerstrokecolor = basicColors[10])
vline!([xBar], color =basicColors[10], linestyle=:dash, label = "")
savefig(figFolder*"NormalMargLike.pdf")

