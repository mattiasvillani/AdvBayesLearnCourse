using Plots
using Distributions
using LaTeXStrings
using ColorSchemes
import ColorSchemes: Paired_12; colors = Paired_12
basicColors = Paired_12[[1,2,3,4,7,8]]

default(legend = true)
figFolder = "/home/mv/Dropbox/Teaching/AdvBayesLearnCourse/Slides/Figures/"

function logmarglike(xBar, n, σ, κ0)
        return -(1/2)*log(κ0/(κ0 + n)) - ((n*xBar^2)/(2*σ^2))*(n/(κ0 + n))
end




κs = 0.0001:0.01:100
BFs = [exp(logmarglike(xBar, n, σ, κ0)) for κ0 in κs]
plot(κs, BFs)

# Data
xBar = 0.8
n = 10
σ = 1
μ0 = 0
κ0 = 0.1
θGrid = -3:0.01:3
margDensH0 = pdf.(Normal(0,σ*√(1/n)), θGrid)
p = plot(xlabel = L"\theta", ylabel = "Marginal likelihood")
plot!(p, θGrid, margDensH0, color = basicColors[6], lw = 2, label = L"H_0")
for (count, κ0) ∈ enumerate([0.1 1 10 100])
        margDensH1 = pdf.(Normal(μ0,σ*√(1/n + 1/κ0)), θGrid)
        plot!(p, θGrid, margDensH1, color = basicColors[count], lw = 2,
                label = LaTeXString("\$\\kappa_0 = $κ0\$"))
end
scatter!([xBar], [0], seriestype = :scatter, color = "red", label = "xBar")
savefig(figFolder*"NormalMargLike.svg")
