## Process results
# Running this script will output svg files to figs&tables. For `BReuss_dh.svg`
#  and `BVoigt_dh.svg` we output to figs&tables/raw and slightly adjust some of
#  the figure features (e.g., fine tune the position of lines to structures etc.)

cd("./res")
using Pkg; Pkg.activate(".")
using PiezoResults
using DelimitedFiles, CairoMakie, ColorSchemes, Colors, LaTeXStrings,
  CategoricalArrays, DataFrames, LaTeXTabulars, JSON3, JLD2, Format

CairoMakie.activate!(type = "svg")
colors = ColorSchemes.seaborn_deep.colors;

## Load data
# Only render and build if need to re-generate images in data/structures/raw AND/OR jld2. Requires raw data.
# render_structures(result_path="/media/zachary/DATA/Data/Piezo/BReuss/")
# render_structures(result_path="/media/zachary/DATA/Data/Piezo/BVoigt/")
# build_max_bh_dataframe(;path=["/media/zachary/DATA/Data/Piezo/BVoigt/","/media/zachary/DATA/Data/Piezo/BReuss/"])
data = load("$(@__DIR__)/../data/max_bh_results.jld2")["data"]

## Plot data
BVoigt_dh = build_bh_dh_plot(data[data.Initial_struc .== "Regular",:],"Bh_Voigt";show_structures=true,structure_size=80,dh_interval=(0.16,0.25))
BReuss_dh = build_bh_dh_plot(data[data.Initial_struc .== "Regular",:],"Bh_Reuss";show_structures=true,structure_size=80,dh_interval=(0.16,0.25))

mkpath((@__DIR__)*"/../figs&tables/raw");
save((@__DIR__)*"/../figs&tables/raw/Figure4_raw.svg",BVoigt_dh)
save((@__DIR__)*"/../figs&tables/raw/Figure5_raw.svg",BReuss_dh)

## Plot all data (regular initial structure)
comb_foms = build_BMeas_FOM_plots(data[data.Initial_struc .== "Regular" .&& data.Req_dh .!= 0.4,:],["dh","gh","kh"])
save((@__DIR__)*"/../figs&tables/raw/Figure6_raw.svg",comb_foms)

## Histories
# This is a convinence operation to place history plots within the data frame.
#   Do not save the resulting data frame, the output will be massive. You are
#   able to save the resulting figures seperately.

## Full dataset
# build_history_plots!(data); # full dataset

## Partial dataset
fontsize = 0.8*30
data_subset = data[data.Initial_struc .== "Regular" .&& data.Req_dh .== 0.29,:]
build_history_plots!(data_subset;show_initial_struc=false,fontsize)
histfig1 = first.(getproperty.((data_subset[data_subset.Objective .== "Bh_Voigt",:],),[:ID,:History_Plot])); histfig1[2]
save((@__DIR__)*"/../figs&tables/Figure3.svg",histfig1[2])

## Base material comparisons
ε_0 = 8.854e-12;
C_base = [12.0400e10  7.52000e10  7.51000e10  0.0        0.0        0.0
      7.52000e10  12.0400e10  7.51000e10  0.0        0.0        0.0
      7.51000e10  7.51000e10  11.0900e10  0.0        0.0        0.0
      0.0         0.0         0.0          2.1000e10  0.0        0.0
      0.0         0.0         0.0          0.0        2.1000e10  0.0
      0.0         0.0         0.0          0.0        0.0       2.30e10]
e_base = [0.0       0.0       0.0        0.0       12.30000  0.0
            0.0       0.0       0.0        12.30000  0.0       0.0
            -5.40000  -5.40000   15.80000   0.0       0.0       0.0]
K_base = [540*ε_0     0        0
              0      540*ε_0    0
              0        0    830*ε_0]

gh_base = PiezoResults.compute_gh(C_base,e_base,K_base)
kh_base = PiezoResults.compute_kh(C_base,e_base,K_base)
S_base = inv(C_base);
d = e_base*S_base;
dh_base = d[3,1]+d[3,2]+d[3,3]
B_Voigt_base = 1/9*(C_base[1,1]+C_base[2,2]+C_base[3,3]+2*(C_base[1,2]+C_base[1,3]+C_base[2,3]))
B_Reuss_base = (1/sum(S_base[i,j] for i = 1:3, j = 1:3))

## max BReuss dh029 vol050
maxbh_dh029_vol050 = data[data.Objective .== "Bh_Reuss" .&& data.Initial_struc .== "Regular" .&& data.Req_dh .== 0.29,:]
gh_maxbh_dh029_vol050 = PiezoResults.compute_gh(maxbh_dh029_vol050.C[1],maxbh_dh029_vol050.e[1],maxbh_dh029_vol050.K[1])
kh_maxbh_dh029_vol050 = PiezoResults.compute_kh(maxbh_dh029_vol050.C[1],maxbh_dh029_vol050.e[1],maxbh_dh029_vol050.K[1])
dh_maxbh_dh029_vol050 = maxbh_dh029_vol050.dh[1]*10^-9
B_Voigt_maxbh_dh029_vol050 = maxbh_dh029_vol050.Bh_Voigt[1]
B_Reuss_maxbh_dh029_vol050 = maxbh_dh029_vol050.Bh_Reuss[1]

dh_maxbh_dh029_vol050/dh_base
gh_maxbh_dh029_vol050/gh_base
kh_maxbh_dh029_vol050/kh_base
B_Voigt_maxbh_dh029_vol050/B_Voigt_base
B_Reuss_maxbh_dh029_vol050/B_Reuss_base

## max BVoigt dh029 vol050
maxbhVoigt_dh029_vol050 = data[data.Objective .== "Bh_Voigt" .&& data.Initial_struc .== "Regular" .&& data.Req_dh .== 0.29,:]
gh_maxbhVoigt_dh029_vol050 = PiezoResults.compute_gh(maxbhVoigt_dh029_vol050.C[1],maxbhVoigt_dh029_vol050.e[1],maxbhVoigt_dh029_vol050.K[1])
kh_maxbhVoigt_dh029_vol050 = PiezoResults.compute_kh(maxbhVoigt_dh029_vol050.C[1],maxbhVoigt_dh029_vol050.e[1],maxbhVoigt_dh029_vol050.K[1])
dh_maxbhVoigt_dh029_vol050 = maxbhVoigt_dh029_vol050.dh[1]*10^-9
B_Voigt_maxbhVoigt_dh029_vol050 = maxbhVoigt_dh029_vol050.Bh_Voigt[1]
B_Reuss_maxbhVoigt_dh029_vol050 = maxbhVoigt_dh029_vol050.Bh_Reuss[1]

dh_maxbhVoigt_dh029_vol050/dh_base
gh_maxbhVoigt_dh029_vol050/gh_base
kh_maxbhVoigt_dh029_vol050/kh_base
B_Voigt_maxbhVoigt_dh029_vol050/B_Voigt_base
B_Reuss_maxbhVoigt_dh029_vol050/B_Reuss_base