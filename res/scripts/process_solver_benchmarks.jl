using Pkg; Pkg.activate(".")
using PiezoResults
using DelimitedFiles, CairoMakie, ColorSchemes, Colors, LaTeXStrings,
  CategoricalArrays, DataFrames, LaTeXTabulars

CairoMakie.activate!(type = "svg")
colors = ColorSchemes.seaborn_deep.colors;

## Load data
# build_benchmark_dataframe() # Run if processing and saving new jld2
data = load("$(@__DIR__)/../data/benchmark_data.jld2")["data"]

schur_data_1eM2 = sort(data[data.Solver .== "Schur" .&& data.rtol .== 1e-2 ,:],:N);
tricg_data_1eM2 = sort(data[data.Solver .== "TriCG" .&& data.rtol .== 1e-2 ,:],:N);
blockGMRES_data_1eM2 = sort(data[data.Solver .== "BlockGMRES" .&& data.rtol .== 1e-2 ,:],:N);
schur_data_1eM4 = sort(data[data.Solver .== "Schur" .&& data.rtol .== 1e-4 ,:],:N);
tricg_data_1eM4 = sort(data[data.Solver .== "TriCG" .&& data.rtol .== 1e-4 ,:],:N);
blockGMRES_data_1eM4 = sort(data[data.Solver .== "BlockGMRES" .&& data.rtol .== 1e-4 ,:],:N);
schur_data_1eM8 = sort(data[data.Solver .== "Schur" .&& data.rtol .== 1e-8 ,:],:N);
tricg_data_1eM8 = sort(data[data.Solver .== "TriCG" .&& data.rtol .== 1e-8 ,:],:N);
blockGMRES_data_1eM8 = sort(data[data.Solver .== "BlockGMRES" .&& data.rtol .== 1e-8 ,:],:N);

## Weak scaling
fig1 = with_theme(theme_latexfonts(), palette=(color=colors,markercolor=color,patchcolor=color)) do
  fig = Figure(fontsize = 20, markersize = 24)
  ax = Axis(fig[2, 1], xlabel = "Number of processors", ylabel="Time (sec)",xticks=schur_data_1eM4.N,
        xtickalign=1,ytickalign=1,xminorticksvisible=true,yminorticksvisible=true,
        xminortickalign=1,yminortickalign=1,xminorticks=IntervalsBetween(5),yminorticks=IntervalsBetween(5))
  scatter!(ax,blockGMRES_data_1eM8.N,blockGMRES_data_1eM8.time,label = L"DBP‑GMRES ($\varepsilon_\mathrm{rel}=10^{-8}$)",color=colors[1],marker=:xcross)
  scatter!(ax,blockGMRES_data_1eM4.N,blockGMRES_data_1eM4.time,label = L"DBP‑GMRES ($\varepsilon_\mathrm{rel}=10^{-4}$)",color=colors[2],marker=:xcross)
  scatter!(ax,blockGMRES_data_1eM2.N,blockGMRES_data_1eM2.time,label = L"DBP‑GMRES ($\varepsilon_\mathrm{rel}=10^{-2}$)",color=colors[7],marker=:xcross)
  scatter!(ax,tricg_data_1eM8.N,tricg_data_1eM8.time,label = L"TriCG ($\varepsilon_\mathrm{rel}=10^{-8}$)",strokecolor=colors[1],strokewidth=3,color=(colors[1],0.0))
  scatter!(ax,tricg_data_1eM4.N,tricg_data_1eM4.time,label = L"TriCG ($\varepsilon_\mathrm{rel}=10^{-4}$)",strokecolor=colors[2],strokewidth=3,color=(colors[2],0.0))
  scatter!(ax,tricg_data_1eM2.N,tricg_data_1eM2.time,label = L"TriCG ($\varepsilon_\mathrm{rel}=10^{-2}$)",strokecolor=colors[7],strokewidth=3,color=(colors[7],0.0))
  scatter!(ax,schur_data_1eM8.N,schur_data_1eM8.time,label = L"SCP‑GMRES ($\varepsilon_\mathrm{rel}=10^{-8}$)",color=colors[1],marker=:star4)
  scatter!(ax,schur_data_1eM4.N,schur_data_1eM4.time,label = L"SCP‑GMRES ($\varepsilon_\mathrm{rel}=10^{-4}$)",color=colors[2],marker=:star4)
  scatter!(ax,schur_data_1eM2.N,schur_data_1eM2.time,label = L"SCP‑GMRES ($\varepsilon_\mathrm{rel}=10^{-2}$)",color=colors[7],marker=:star4)
  Legend(fig[1,1],ax,orientation = :horizontal,nbanks=3, tellwidth = true)
  resize_to_layout!(fig)
  fig
end

save((@__DIR__)*"/../figs&tables/Figure2.svg",fig1)

## Table
_N = schur_data_1eM4.N
table_data = [
  Rule(:h),
  ["Method","",MultiColumn(length(_N), :l, "Number of processors")], # \cline{3-7}
  ["","",_N...],
  Rule(:h),
  ["DBP-GMRES",L"$\boldsymbol{v}\mapsto A\backslash \boldsymbol{v}$ Iters.",blockGMRES_data_1eM8.U_its...],
  [L"($\varepsilon_\mathrm{rel}=10^{-8}$)",L"$\boldsymbol{u}\mapsto -C\backslash \boldsymbol{u}$ Iters.",blockGMRES_data_1eM8.phi_its...],
  ["","Global Iters.",blockGMRES_data_1eM8.outer_its...],
  ["DBP-GMRES",L"$\boldsymbol{v}\mapsto A\backslash \boldsymbol{v}$ Iters.",blockGMRES_data_1eM4.U_its...],
  [L"($\varepsilon_\mathrm{rel}=10^{-4}$)",L"$\boldsymbol{u}\mapsto -C\backslash \boldsymbol{u}$ Iters.",blockGMRES_data_1eM4.phi_its...],
  ["","Global Iters.",blockGMRES_data_1eM4.outer_its...],
  ["DBP-GMRES",L"$\boldsymbol{v}\mapsto A\backslash \boldsymbol{v}$ Iters.",blockGMRES_data_1eM2.U_its...],
  [L"($\varepsilon_\mathrm{rel}=10^{-2}$)",L"$\boldsymbol{u}\mapsto -C\backslash \boldsymbol{u}$ Iters.",blockGMRES_data_1eM2.phi_its...],
  ["","Global Iters.",blockGMRES_data_1eM2.outer_its...],
  ["TriCG",L"$\boldsymbol{v}\mapsto A\backslash \boldsymbol{v}$ Iters.",tricg_data_1eM8.U_its...],
  [L"($\varepsilon_\mathrm{rel}=10^{-8}$)",L"$\boldsymbol{u}\mapsto -C\backslash \boldsymbol{u}$ Iters.",tricg_data_1eM8.phi_its...],
  ["","Global Iters.",tricg_data_1eM8.outer_its...],
  ["TriCG",L"$\boldsymbol{v}\mapsto A\backslash \boldsymbol{v}$ Iters.",tricg_data_1eM4.U_its...],
  [L"($\varepsilon_\mathrm{rel}=10^{-4}$)",L"$\boldsymbol{u}\mapsto -C\backslash \boldsymbol{u}$ Iters.",tricg_data_1eM4.phi_its...],
  ["","Global Iters.",tricg_data_1eM4.outer_its...],
  ["TriCG",L"$\boldsymbol{v}\mapsto A\backslash \boldsymbol{v}$ Iters.",tricg_data_1eM2.U_its...],
  [L"($\varepsilon_\mathrm{rel}=10^{-2}$)",L"$\boldsymbol{u}\mapsto -C\backslash \boldsymbol{u}$ Iters.",tricg_data_1eM2.phi_its...],
  ["","Global Iters.",tricg_data_1eM2.outer_its...],
  ["SCP-GMRES",L"$\boldsymbol{v}\mapsto A\backslash \boldsymbol{v}$ Iters.",schur_data_1eM8.U_its...],
  [L"($\varepsilon_\mathrm{rel}=10^{-8}$)",L"$\boldsymbol{u}\mapsto \bar{S}\backslash \boldsymbol{u}$ Iters.",schur_data_1eM8.phi_its...],
  ["","Global Iters.",schur_data_1eM8.outer_its...],
  ["SCP-GMRES",L"$\boldsymbol{v}\mapsto A\backslash \boldsymbol{v}$ Iters.",schur_data_1eM4.U_its...],
  [L"($\varepsilon_\mathrm{rel}=10^{-4}$)",L"$\boldsymbol{u}\mapsto \bar{S}\backslash \boldsymbol{u}$ Iters.",schur_data_1eM4.phi_its...],
  ["","Global Iters.",schur_data_1eM4.outer_its...],
  ["SCP-GMRES",L"$\boldsymbol{v}\mapsto A\backslash \boldsymbol{v}$ Iters.",schur_data_1eM2.U_its...],
  [L"($\varepsilon_\mathrm{rel}=10^{-2}$)",L"$\boldsymbol{u}\mapsto \bar{S}\backslash \boldsymbol{u}$ Iters.",schur_data_1eM2.phi_its...],
  ["","Global Iters.",schur_data_1eM2.outer_its...],
]

# Note, need to add  \cline{3-7} after first instance of `\\`
# also add this after each method (roughly every third instance of `\\`)
# Add \hline after last instance of `\\`
latex_tabular((@__DIR__)*"/../figs&tables/Table1.tex",Tabular('l'^(length(_N)+2)),table_data)