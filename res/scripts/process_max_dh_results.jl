cd("./res")
using Pkg; Pkg.activate(".")
using PiezoResults
using DelimitedFiles, CairoMakie, ColorSchemes, Colors, LaTeXStrings,
  CategoricalArrays, DataFrames, LaTeXTabulars, JSON3, JLD2, Format

CairoMakie.activate!(type = "svg")
colors = ColorSchemes.seaborn_deep.colors

## Load data
# Only render and build if need to re-generate images in data/structures/raw AND/OR jld2. Requires raw data.
# render_structures(result_path="/media/zachary/DATA/Data/Piezo/dh/")
# build_max_dh_dataframe(;path=["/media/zachary/DATA/Data/Piezo/dh/"])
data = load("$(@__DIR__)/../data/max_dh_results.jld2")["data"]
d1 = data[3:3,:]; d2 = data[2:2,:]; d3 = data[1:1,:];
font_size_w_init_struc = 30

_tmp = [d1.History[1].J;d2.History[1].J;d3.History[1].J]
_ymin = -maximum(_tmp); _ymin *= 0.99
_ymax = -minimum(_tmp); _ymax *= 1.01
ylims = (:manual,_ymin,_ymax)

build_history_plots_dh!(d1;ylims,show_initial_struc=false,show_labels=true,left_label_only = true, skip_i_letters=0,fontsize=font_size_w_init_struc*0.8)
build_history_plots_dh!(d2;ylims,show_initial_struc=false,show_labels=true,left_label_only = true, skip_i_letters=1,fontsize=font_size_w_init_struc*0.8)
build_history_plots_dh!(d3;ylims,show_initial_struc=false,show_labels=true,left_label_only = true, skip_i_letters=2,fontsize=font_size_w_init_struc*0.8)

fig1 = d1.History_Plot[1]
fig2 = d2.History_Plot[1]
fig3 = d3.History_Plot[1]

save((@__DIR__)*"/../figs&tables/Figure7a.svg",fig1)
save((@__DIR__)*"/../figs&tables/Figure7b.svg",fig2)
save((@__DIR__)*"/../figs&tables/Figure7c.svg",fig3)

## Table of properties
table_data = [
  ["Fig.",L"Required $\bar{C}_{zzzz}$ (Pa)",L"$\max\bar{d}_h$ (C/N)",L"\mathrm{Vol}",L"$\bar{B}_\mathrm{Voigt}$ (Pa)",L"$\bar{B}_\mathrm{Reuss}$ (Pa)",L"$\bar{g}_h$ (Vm/N)",L"\bar{k}_h"],
  Rule(:h),
  ["(a)";round.([d1.Req_C33[1]*10^10,d1.dh[1],d1.Vol[1],d1.Bh_Voigt[1],d1.Bh_Reuss[1],PiezoResults.compute_gh(d1.C[1],d1.e[1],d1.K[1]),PiezoResults.compute_kh(d1.C[1],d1.e[1],d1.K[1])];sigdigits=3)],
  Rule(:h),
  ["(b)";round.([d2.Req_C33[1]*10^10,d2.dh[1],d2.Vol[1],d2.Bh_Voigt[1],d2.Bh_Reuss[1],PiezoResults.compute_gh(d2.C[1],d2.e[1],d2.K[1]),PiezoResults.compute_kh(d2.C[1],d2.e[1],d2.K[1])];sigdigits=3)],
  Rule(:h),
  ["(c)";round.([d3.Req_C33[1]*10^10,d3.dh[1],d3.Vol[1],d3.Bh_Voigt[1],d3.Bh_Reuss[1],PiezoResults.compute_gh(d3.C[1],d3.e[1],d3.K[1]),PiezoResults.compute_kh(d3.C[1],d3.e[1],d3.K[1])];sigdigits=3)],
  Rule(:h)
]

latex_tabular((@__DIR__)*"/../figs&tables/Table2.tex",Tabular("l|"*'l'^(7)),table_data)