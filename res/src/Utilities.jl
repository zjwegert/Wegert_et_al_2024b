alpharize(a, colors) = RGBAf.(colors, a)

function render_structures(;result_path,output_folder="$(@__DIR__)/../data/structures/raw",pvpython_cmd="pvpython")
  if Sys.iswindows()
    run(`cmd /c $pvpython_cmd $(@__DIR__)/render_structures.py $result_path $output_folder`)
  else
    run(`$pvpython_cmd $(@__DIR__)/render_structures.py $result_path $output_folder`)
  end
end

# NOTE: We multiple dh by 10 as the output is in normalised units which has scalign 10^-9.
#  Multiplying by 10 allows ticks to be order 1 and label says 10^-10
function build_bh_dh_plot(subdata,objective;
    show_structures=true,
    structure_size=100,
    structure_path="$(@__DIR__)/../data/structures/raw/",
    dh_interval = (0.14,0.29),
    colors=ColorSchemes.seaborn_deep.colors,
    fontsize = 18,
    markersize = 24
    )
  init_strucs = unique(subdata.Initial_struc); @assert length(init_strucs)==1 "Expected data pertaining to the same initial structure"
  @assert objective ∈ ["Bh_Reuss","Bh_Voigt"] "Expected objective to be either Bh_Reuss or Bh_Voigt"
  opt_data = subdata[subdata.Objective .== objective,:]
  aux_data = subdata[subdata.Objective .!= objective,:]
  measure_name = split(objective,"Bh_")[2];
  aux_measure_name = split(first(aux_data.Objective),"Bh_")[2]
  figure = with_theme(theme_latexfonts(), palette=(color=colors,markercolor=color,patchcolor=color)) do
    fig = Figure(fontsize = fontsize, markersize = markersize)
    ax = Axis(fig[1, 1], xlabel = L"$\bar{B}_{\mathrm{%$measure_name}}$ ($\times 10^{10}$ N/m^2)", ylabel=L"$\bar{d}_h$ ($\times 10^{-10}$ C/N)",
          xtickalign=1,ytickalign=1,xminorticksvisible=true,yminorticksvisible=true,
          xminortickalign=1,yminortickalign=1,xminorticks=IntervalsBetween(5),yminorticks=IntervalsBetween(5))
    if objective == "Bh_Voigt"
      scatter!(ax,getproperty(opt_data,Symbol(first(opt_data.Objective)))/10^10,opt_data.dh*10,label = L"$\max_\Omega~\bar{B}_{\mathrm{%$measure_name}}$",strokecolor=colors[1],strokewidth=3,color=(colors[1],0.0))
      scatter!(ax,getproperty(aux_data,Symbol(first(opt_data.Objective)))/10^10,aux_data.dh*10,label = L"$\max_\Omega~\bar{B}_{\mathrm{%$aux_measure_name}}$",marker=:xcross,color=colors[2])
    else
      scatter!(ax,getproperty(aux_data,Symbol(first(opt_data.Objective)))/10^10,aux_data.dh*10,label = L"$\max_\Omega~\bar{B}_{\mathrm{%$aux_measure_name}}$",strokecolor=colors[1],strokewidth=3,color=(colors[1],0.0))
      scatter!(ax,getproperty(opt_data,Symbol(first(opt_data.Objective)))/10^10,opt_data.dh*10,label = L"$\max_\Omega~\bar{B}_{\mathrm{%$measure_name}}$",marker=:xcross,color=colors[2])
    end

    if show_structures
      g_axis_data = Dict()
      ga = fig[0,1:2] = GridLayout() # Change 1:2 to 1 to remove overlap
      for (i,dt) in enumerate(eachrow(sort!(opt_data[opt_data.Req_dh .>= dh_interval[2],:],[:Req_dh],rev=true)))
        ax2 = Axis(ga[1,i],aspect = DataAspect(),yticksvisible=false,xticksvisible=false, width = structure_size, height = structure_size)
        img_SOLID = load(joinpath(structure_path,dt.ID*"_STRUCTURE.png"));
        image!(ax2, rotr90(img_SOLID))
        hidedecorations!(ax2)
        hidespines!(ax2)
        g_axis_data[dt] = ax2
      end
      colgap!(ga, 0)
      gb = fig[1,2] = GridLayout()
      for (i,dt) in enumerate(eachrow(sort!(opt_data[dh_interval[1] .< opt_data.Req_dh .< dh_interval[2],:],[:Req_dh],rev=true)))
        ax2 = Axis(gb[i,1],aspect = DataAspect(),yticksvisible=false,xticksvisible=false, width = structure_size, height = structure_size)
        img_SOLID = load(joinpath(structure_path,dt.ID*"_STRUCTURE.png"));
        image!(ax2, rotr90(img_SOLID))
        hidedecorations!(ax2)
        hidespines!(ax2)
        g_axis_data[dt] = ax2
      end
      rowgap!(gb, 0)
      gc = fig[2,1:2] = GridLayout() # Change 1:2 to 1 to remove overlap
      for (i,dt) in enumerate(eachrow(sort!(opt_data[opt_data.Req_dh .<= dh_interval[1],:],[:Req_dh],rev=false)))
        ax2 = Axis(gc[1,i],aspect = DataAspect(),yticksvisible=false,xticksvisible=false, width = structure_size, height = structure_size)
        img_SOLID = load(joinpath(structure_path,dt.ID*"_STRUCTURE.png"));
        image!(ax2, rotr90(img_SOLID))
        hidedecorations!(ax2)
        hidespines!(ax2)
        g_axis_data[dt] = ax2
      end
      colgap!(gc, 0)
    end

    # Legend(fig[3-(1-show_structures),1],ax,orientation = :horizontal)
    axislegend(ax; position = :lb)
    resize_to_layout!(fig)

    if show_structures
      x = Float64[]; y = Float64[]
      u = Float64[]; v = Float64[]
      for dt in keys(g_axis_data)
        _xy = Makie.shift_project(ax.scene,fig.scene,(getproperty(dt,Symbol(dt.Objective))*10^-10,dt.Req_dh*10))
        _uv = Makie.shift_project(g_axis_data[dt].scene,fig.scene,(800,800)) - _xy
        push!(x,_xy[1]); push!(y,_xy[2])
        push!(u,_uv[1]); push!(v,_uv[2])
      end
      arrows!(fig.scene,x,y,u,v,space=:pixel,color=(:grey, 0.5))
    end

    fig
  end
  return figure
end

function build_bh_dh_plot_alternate(subdata,objective;
    show_structures=true,
    structure_size=100,
    structure_path="$(@__DIR__)/../data/structures/raw/",
    dh_interval = (0.14,0.29),
    colors=ColorSchemes.seaborn_deep.colors,
    fontsize = 18,
    markersize = 24
    )
  xtickmin=Base._round_step(minimum([subdata.Bh_Reuss/10^10;subdata.Bh_Voigt/10^10]),0.5,RoundDown)-0.1;
  xtickmax=Base._round_step(maximum([subdata.Bh_Reuss/10^10;subdata.Bh_Voigt/10^10]),0.5,RoundUp);
  init_strucs = unique(subdata.Initial_struc); @assert length(init_strucs)==1 "Expected data pertaining to the same initial structure"
  @assert objective ∈ ["Bh_Reuss","Bh_Voigt"] "Expected objective to be either Bh_Reuss or Bh_Voigt"
  opt_data = subdata[subdata.Objective .== objective,:]
  function kwargs(col)
    if objective == "Bh_Voigt"
      return (strokecolor=col,strokewidth=3,color=(col,0.0))
    else
      return (marker=:xcross,color=col)
    end
  end
  figure = with_theme(theme_latexfonts(), palette=(color=colors,markercolor=color,patchcolor=color)) do
    fig = Figure(fontsize = fontsize, markersize = markersize)
    ax = Axis(fig[1, 1], xlabel = L"$\bar{B}_{\mathrm{Voigt}}$ ($\times 10^{10}$ N/m^2)", ylabel=L"$\bar{d}_h$ ($\times 10^{-10}$ C/N)",
          xticklabelcolor = colors[1], xlabelcolor = colors[1],xtickalign=1,ytickalign=1,xminorticksvisible=true,yminorticksvisible=true,
          xminortickalign=1,yminortickalign=1,xminorticks=IntervalsBetween(5),yminorticks=IntervalsBetween(5),
          xticks = collect(0:0.5:xtickmax))
    ax_2 = Axis(fig[1, 1], xlabel = L"$\bar{B}_{\mathrm{Reuss}}$ ($\times 10^{10}$ N/m^2)", xlabelcolor = colors[2],
          xaxisposition=:top,xticklabelcolor = colors[2],xtickalign=1,ytickalign=1,xminorticksvisible=true,yminorticksvisible=true,
          xminortickalign=1,yminortickalign=1,xminorticks=IntervalsBetween(5),yminorticks=IntervalsBetween(5),
          xticks = collect(0:0.5:xtickmax))
    xlims!(ax,xtickmin,xtickmax)
    xlims!(ax_2,xtickmin,xtickmax)
    hidespines!(ax_2)
    hideydecorations!(ax_2)
    scatter!(ax,opt_data.Bh_Voigt/10^10,opt_data.dh*10;kwargs(colors[1])...)
    scatter!(ax_2,opt_data.Bh_Reuss/10^10,opt_data.dh*10;kwargs(colors[2])...)

    if show_structures
      g_axis_data = Dict()
      ga = fig[0,1:2] = GridLayout() # Change 1:2 to 1 to remove overlap
      for (i,dt) in enumerate(eachrow(sort!(opt_data[opt_data.Req_dh .>= dh_interval[2],:],[:Req_dh],rev=true)))
        ax2 = Axis(ga[1,i],aspect = DataAspect(),yticksvisible=false,xticksvisible=false, width = structure_size, height = structure_size)
        img_SOLID = load(joinpath(structure_path,dt.ID*"_STRUCTURE.png"));
        image!(ax2, rotr90(img_SOLID))
        hidedecorations!(ax2)
        hidespines!(ax2)
        g_axis_data[dt] = ax2
      end
      colgap!(ga, 0)
      gb = fig[1,2] = GridLayout()
      for (i,dt) in enumerate(eachrow(sort!(opt_data[dh_interval[1] .< opt_data.Req_dh .< dh_interval[2],:],[:Req_dh],rev=true)))
        ax2 = Axis(gb[i,1],aspect = DataAspect(),yticksvisible=false,xticksvisible=false, width = structure_size, height = structure_size)
        img_SOLID = load(joinpath(structure_path,dt.ID*"_STRUCTURE.png"));
        image!(ax2, rotr90(img_SOLID))
        hidedecorations!(ax2)
        hidespines!(ax2)
        g_axis_data[dt] = ax2
      end
      rowgap!(gb, 0)
      gc = fig[2,1:2] = GridLayout() # Change 1:2 to 1 to remove overlap
      for (i,dt) in enumerate(eachrow(sort!(opt_data[opt_data.Req_dh .<= dh_interval[1],:],[:Req_dh],rev=false)))
        ax2 = Axis(gc[1,i],aspect = DataAspect(),yticksvisible=false,xticksvisible=false, width = structure_size, height = structure_size)
        img_SOLID = load(joinpath(structure_path,dt.ID*"_STRUCTURE.png"));
        image!(ax2, rotr90(img_SOLID))
        hidedecorations!(ax2)
        hidespines!(ax2)
        g_axis_data[dt] = ax2
      end
      colgap!(gc, 0)
    end

    resize_to_layout!(fig)

    if show_structures
      x = Float64[]; y = Float64[]
      u = Float64[]; v = Float64[]
      for dt in keys(g_axis_data)
        _xy = Makie.shift_project(ax.scene,fig.scene,(getproperty(dt,Symbol(dt.Objective))*10^-10,dt.Req_dh*10))
        _uv = Makie.shift_project(g_axis_data[dt].scene,fig.scene,(800,800)) - _xy
        push!(x,_xy[1]); push!(y,_xy[2])
        push!(u,_uv[1]); push!(v,_uv[2])
      end
      arrows!(fig.scene,x,y,u,v,space=:pixel,color=(:grey, 0.5))
    end

    fig
  end
  return figure
end

function compute_gh(C,e,K)
  d = e*inv(C);
  g = inv(K+d*C*d')*d;
  return g[3,1]+g[3,2]+g[3,3];
end

function compute_kh(C,e,K)
  d = e*inv(C);
  Ksigma = K+d*C*d';
  dh = d[3,1]+d[3,2]+d[3,3];
  S = inv(C);
  Sh = sum(S[i,j] for i = 1:3, j = 1:3)
  return sqrt(dh^2/(Ksigma[3,3]*Sh));
end

function build_BMeas_FOM_plots(data,fom::Vector{String};
    colors=ColorSchemes.seaborn_deep.colors,
    fontsize = 18,
    markersize = 24,
    aspect = 2)
  @assert all(x-> x ∈ ["dh","gh","kh"],fom)
  init_strucs = unique(data.Initial_struc); @assert length(init_strucs)==1 "Expected data pertaining to the same initial structure"
  data_reuss = data[data.Objective .== "Bh_Reuss",:]
  data_voigt = data[data.Objective .== "Bh_Voigt",:]
  count = 0
  letters = collect('a':'z')
  axs = []
  _fig = with_theme(theme_latexfonts(), palette=(color=colors,markercolor=color,patchcolor=color)) do
    fig = Figure(fontsize = fontsize, markersize = markersize)
    xtickmin=0.0#Base._round_step(minimum([data.Bh_Reuss/10^10;data.Bh_Voigt/10^10]),0.5,RoundDown)-0.1;
    xtickmax=Base._round_step(maximum([data.Bh_Reuss/10^10;data.Bh_Voigt/10^10]),0.5,RoundUp);
    if "dh" ∈ fom
      count += 1
      show_xlabel_and_ticks = count > 1 ? false : true
      _ylabel = L"$\bar{d}_h$ ($\times 10^{-10}$ C/N)"
      _yscale = identity
      y_data_voigt = data_voigt.dh*10 # We do this because scaling is 10^-9 on data, so mult by 10 to get ticks order 1 and label 10^-10
      y_data_reuss = data_reuss.dh*10
      ax = Axis(fig[count, 1], xlabel = L"$\bar{B}_{\mathrm{Reuss}}$ ($\times 10^{10}$ N/m^2)", ylabel=_ylabel,
            xticklabelcolor = (colors[1],show_xlabel_and_ticks ? 1 : 0), xlabelcolor = (colors[1],show_xlabel_and_ticks ? 1 : 0),
            xtickalign=1,ytickalign=1,xminorticksvisible=true,yminorticksvisible=true,
            xminortickalign=1,yminortickalign=1,xminorticks=IntervalsBetween(5),yminorticks=IntervalsBetween(5),
            xticks = collect(0:0.5:xtickmax),yscale=_yscale,height=200)
      ax2 = Axis(fig[count, 1], xlabel = L"$\bar{B}_{\mathrm{Voigt}}$ ($\times 10^{10}$ N/m^2)",
            xlabelcolor = (colors[2],show_xlabel_and_ticks ? 1 : 0),
            xaxisposition=:top,xticklabelcolor = (colors[2],show_xlabel_and_ticks ? 1 : 0),
            xtickalign=1,ytickalign=1,xminorticksvisible=true,yminorticksvisible=true,
            xminortickalign=1,yminortickalign=1,xminorticks=IntervalsBetween(5),yminorticks=IntervalsBetween(5),
            xticks = collect(0:0.5:xtickmax),yscale=_yscale,height=200)
      push!(axs,ax,ax2)
      xlims!(ax,xtickmin,xtickmax)
      xlims!(ax2,xtickmin,xtickmax)
      hidespines!(ax2)
      hideydecorations!(ax2)
      scatter!(ax2,data_voigt.Bh_Voigt/10^10,y_data_voigt,strokecolor=colors[2],strokewidth=3,color=(colors[2],0.0))
      scatter!(ax,data_voigt.Bh_Reuss/10^10,y_data_voigt,strokecolor=colors[1],strokewidth=3,color=(colors[1],0.0))
      scatter!(ax2,data_reuss.Bh_Voigt/10^10,y_data_reuss,marker=:xcross,color=colors[2])
      scatter!(ax,data_reuss.Bh_Reuss/10^10,y_data_reuss,marker=:xcross,color=colors[1])
      axislegend(ax,[[MarkerElement(color = (:black,0.0), marker = :circle, markersize = 24,strokecolor = :black,strokewidth=3)],
        [MarkerElement(color = :black, marker = :xcross, markersize = 24,strokecolor = :black),]],
        [L"$\max_\Omega~\bar{B}_{\mathrm{Voigt}}$",
        L"$\max_\Omega~\bar{B}_{\mathrm{Reuss}}$"],position=:lb)
      length(fom) > 1 ? Label(fig[count, 1, TopLeft()], rich("($(letters[count]))")) : nothing;
    end
    if "gh" ∈ fom
      count += 1
      show_xlabel_and_ticks = count > 1 ? false : true
      _ylabel = L"$\bar{g}_h$ (Vm/N)"
      _yscale = log10
      y_data_voigt = compute_gh.(data_voigt.C,data_voigt.e,data_voigt.K)
      y_data_reuss = compute_gh.(data_reuss.C,data_reuss.e,data_reuss.K)
            ax = Axis(fig[count, 1], xlabel = L"$\bar{B}_{\mathrm{Reuss}}$ ($\times 10^{10}$ N/m^2)", ylabel=_ylabel,
            xticklabelcolor = (colors[1],show_xlabel_and_ticks ? 1 : 0), xlabelcolor = (colors[1],show_xlabel_and_ticks ? 1 : 0),
            xtickalign=1,ytickalign=1,xminorticksvisible=true,yminorticksvisible=true,
            xminortickalign=1,yminortickalign=1,xminorticks=IntervalsBetween(5),yminorticks=IntervalsBetween(5),
            xticks = collect(0:0.5:xtickmax),yscale=_yscale,height=200)
      ax2 = Axis(fig[count, 1], xlabel = L"$\bar{B}_{\mathrm{Voigt}}$ ($\times 10^{10}$ N/m^2)",
            xlabelcolor = (colors[2],show_xlabel_and_ticks ? 1 : 0),
            xaxisposition=:top,xticklabelcolor = (colors[2],show_xlabel_and_ticks ? 1 : 0),
            xtickalign=1,ytickalign=1,xminorticksvisible=true,yminorticksvisible=true,
            xminortickalign=1,yminortickalign=1,xminorticks=IntervalsBetween(5),yminorticks=IntervalsBetween(5),
            xticks = collect(0:0.5:xtickmax),yscale=_yscale,height=200)
      push!(axs,ax,ax2)
      xlims!(ax,xtickmin,xtickmax)
      xlims!(ax2,xtickmin,xtickmax)
      hidespines!(ax2)
      hideydecorations!(ax2)
      scatter!(ax2,data_voigt.Bh_Voigt/10^10,y_data_voigt,strokecolor=colors[2],strokewidth=3,color=(colors[2],0.0))
      scatter!(ax,data_voigt.Bh_Reuss/10^10,y_data_voigt,strokecolor=colors[1],strokewidth=3,color=(colors[1],0.0))
      scatter!(ax2,data_reuss.Bh_Voigt/10^10,y_data_reuss,marker=:xcross,color=colors[2])
      scatter!(ax,data_reuss.Bh_Reuss/10^10,y_data_reuss,marker=:xcross,color=colors[1])
      axislegend(ax,[[MarkerElement(color = (:black,0.0), marker = :circle, markersize = 24,strokecolor = :black,strokewidth=3)],
        [MarkerElement(color = :black, marker = :xcross, markersize = 24,strokecolor = :black),]],
        [L"$\max_\Omega~\bar{B}_{\mathrm{Voigt}}$",
        L"$\max_\Omega~\bar{B}_{\mathrm{Reuss}}$"],position=:lb)
      length(fom) > 1 ? Label(fig[count, 1, TopLeft()], rich("($(letters[count]))")) : nothing;
    end
    if "kh" ∈ fom
      count += 1
      show_xlabel_and_ticks = count > 1 ? false : true
      _ylabel = L"$\bar{k}_h$"
      _yscale = identity
      y_data_voigt = compute_kh.(data_voigt.C,data_voigt.e,data_voigt.K)
      y_data_reuss = compute_kh.(data_reuss.C,data_reuss.e,data_reuss.K)
            ax = Axis(fig[count, 1], xlabel = L"$\bar{B}_{\mathrm{Reuss}}$ ($\times 10^{10}$ N/m^2)", ylabel=_ylabel,
            xticklabelcolor = (colors[1],show_xlabel_and_ticks ? 1 : 0), xlabelcolor = (colors[1],show_xlabel_and_ticks ? 1 : 0),
            xtickalign=1,ytickalign=1,xminorticksvisible=true,yminorticksvisible=true,
            xminortickalign=1,yminortickalign=1,xminorticks=IntervalsBetween(5),yminorticks=IntervalsBetween(5),
            xticks = collect(0:0.5:xtickmax),yscale=_yscale,height=200)
      ax2 = Axis(fig[count, 1], xlabel = L"$\bar{B}_{\mathrm{Voigt}}$ ($\times 10^{10}$ N/m^2)",
            xlabelcolor = (colors[2],show_xlabel_and_ticks ? 1 : 0),
            xaxisposition=:top,xticklabelcolor = (colors[2],show_xlabel_and_ticks ? 1 : 0),
            xtickalign=1,ytickalign=1,xminorticksvisible=true,yminorticksvisible=true,
            xminortickalign=1,yminortickalign=1,xminorticks=IntervalsBetween(5),yminorticks=IntervalsBetween(5),
            xticks = collect(0:0.5:xtickmax),yscale=_yscale,height=200)
      push!(axs,ax,ax2)
      xlims!(ax,xtickmin,xtickmax)
      xlims!(ax2,xtickmin,xtickmax)
      hidespines!(ax2)
      hideydecorations!(ax2)
      scatter!(ax2,data_voigt.Bh_Voigt/10^10,y_data_voigt,strokecolor=colors[2],strokewidth=3,color=(colors[2],0.0))
      scatter!(ax,data_voigt.Bh_Reuss/10^10,y_data_voigt,strokecolor=colors[1],strokewidth=3,color=(colors[1],0.0))
      scatter!(ax2,data_reuss.Bh_Voigt/10^10,y_data_reuss,marker=:xcross,color=colors[2])
      scatter!(ax,data_reuss.Bh_Reuss/10^10,y_data_reuss,marker=:xcross,color=colors[1])
      axislegend(ax,[[MarkerElement(color = (:black,0.0), marker = :circle, markersize = 24,strokecolor = :black,strokewidth=3)],
        [MarkerElement(color = :black, marker = :xcross, markersize = 24,strokecolor = :black),]],
        [L"$\max_\Omega~\bar{B}_{\mathrm{Voigt}}$",
        L"$\max_\Omega~\bar{B}_{\mathrm{Reuss}}$"],position=:lb)
      length(fom) > 1 ? Label(fig[count, 1, TopLeft()], rich("($(letters[count]))")) : nothing;
    end
    resize_to_layout!(fig)
    fig
  end

  yspace = maximum(tight_yticklabel_spacing!, axs)
  xspace = maximum(tight_xticklabel_spacing!, axs)
  for ax in axs
    ax.yticklabelspace = yspace
    ax.xticklabelspace = xspace
  end

  return _fig
end

function build_history_plots!(data;
    structure_path="$(@__DIR__)/../data/structures/raw/",
    colors=ColorSchemes.seaborn_deep.colors,
    new_column_name = "History_Plot",
    show_initial_struc = true,
    show_labels = true,
    left_label_only = false,
    skip_i_letters=0,
    fontsize=20,
    structure_size = 280,
    plots_size = 250,
    nflts = 1,
    ylims::Tuple = (:auto,nothing,nothing) # use (:manual,v1,v2) alternatively
  )
  IDs = data.ID
  ID_to_fig = Dict{String,Figure}()
  letters = collect('a':'z')[1+skip_i_letters:4+skip_i_letters]
  _bold(x) = rich(x)#,font=:bold)
  for ID in IDs
    loc_data = data[data.ID .== ID,:]
    @assert size(loc_data,1) == 1 "We require unique idenitfiers (IDs) for each optimisation solution"
    img_SOLID = load(joinpath(structure_path,first(loc_data.ID)*"_STRUCTURE.png"));
    img_VOID = load(joinpath(structure_path,first(loc_data.ID)*"_STRUCTURE_VOID.png"));
    with_theme(theme_latexfonts(), palette=(color=colors,markercolor=color,patchcolor=color)) do
      f = Figure(fontsize=fontsize, markersize = 24)
      if show_initial_struc
        axes = Axis(f[1,0],title = "Initial stucture",titlecolor=(:black,0.0),aspect = DataAspect(),yticksvisible=false,xticksvisible=false, width = structure_size, height = structure_size)
        Label(f[1, 0, TopLeft()], _bold("($(letters[1]))"),color=(:black,show_labels ? 1.0 : 0.0));
        img_init = load(joinpath(structure_path,loc_data.Initial_struc[1]*"_STRUCTURE.png"));
        image!(axes, rotr90(img_init))
        hidedecorations!(axes)
        hidespines!(axes)
      end
      axes = Axis(f[1,1],title = "Final Solid",titlecolor=(:black,0.0),aspect = DataAspect(),yticksvisible=false,xticksvisible=false, width = structure_size, height = structure_size)
      Label(f[1, 1, TopLeft()], (show_initial_struc ? _bold("($(letters[2]))") : _bold("($(letters[1]))")),color=(:black,show_labels && (~left_label_only || ~show_initial_struc) ? 1.0 : 0.0));
      image!(axes, rotr90(img_SOLID))
      hidedecorations!(axes)
      hidespines!(axes)
      axes = Axis(f[1,2],title = "Final void",titlecolor=(:black,0.0),aspect = DataAspect(),yticksvisible=false,xticksvisible=false, width = structure_size, height = structure_size)
      Label(f[1, 2, TopLeft()], (show_initial_struc ? _bold("($(letters[3]))") : _bold("($(letters[2]))")),color=(:black,show_labels && ~left_label_only ? 1.0 : 0.0));
      image!(axes, rotr90(img_VOID))
      hidedecorations!(axes)
      hidespines!(axes)
      ax1 = Axis(f[1, 3], yticklabelcolor = colors[1], ylabelcolor = colors[1], xlabel = L"\text{Iteration}", ylabel=L"$\bar{B}_{\mathrm{%$(loc_data.Objective[1][4:end])}}(\Omega)\times10^{-10}$",
        xtickalign=1,ytickalign=1,xminorticksvisible=true,yminorticksvisible=true,xgridvisible=true,ygridvisible=true,
        xminortickalign=1,yminortickalign=1,xminorticks=IntervalsBetween(5),yminorticks=IntervalsBetween(5),ytickformat = "{:.2f}", width = plots_size, height = plots_size)
      if ylims[1] == :manual
        ylims!(ax1,ylims[2],ylims[3])
      end
      ax2 = Axis(f[1, 3], yticklabelcolor = colors[2], ylabelcolor = colors[2], yaxisposition = :right, ylabel=L"$C_i(\Omega)$",
        xtickalign=1,ytickalign=1,xminorticksvisible=true,yminorticksvisible=true,xgridvisible=true,ygridvisible=false,
        xminortickalign=1,yminortickalign=1,xminorticks=IntervalsBetween(5),yminorticks=IntervalsBetween(5),ytickformat = "{:.2f}", width = plots_size, height = plots_size,)
      hidespines!(ax2)
      hidexdecorations!(ax2)
      lines!(ax1,-first(loc_data.History).J,color=colors[1],linewidth = 3)
      req_vol = cfmt("%.$(nflts)f",loc_data.Req_Vol[1]);
      req_dh = cfmt("%.$(nflts)f",loc_data.Req_dh[1]*10); # Fix units
      lines!(ax2,first(loc_data.History).Vol,linestyle=:dot,color=colors[2],linewidth = 3,label=L"""\mathrm{Vol}(\Omega)-%$(req_vol)""")
      lines!(ax2,first(loc_data.History).dh,linestyle=:dash,color=colors[2],linewidth = 3,label=L"10^{10}\bar{d}_h(\Omega)-%$(req_dh)")
      Label(f[1, 2, TopRight()], show_initial_struc ? _bold("($(letters[4]))") : _bold("($(letters[3]))"),color=(:black,show_labels && ~left_label_only ? 1.0 : 0.0))#,padding=(-100,0,0,0))
      Legend(f[0,3],ax2,orientation = :horizontal,margin=(0,0,-70,0))
      resize_to_layout!(f)
      ID_to_fig[ID]=f
    end
  end
  data[:,new_column_name] = getindex.((ID_to_fig,),data.ID);
  return data
end

function build_history_plots_dh!(data;
    structure_path="$(@__DIR__)/../data/structures/raw/",
    colors=ColorSchemes.seaborn_deep.colors,
    new_column_name = "History_Plot",
    show_initial_struc = true,
    show_labels = true,
    left_label_only = false,
    skip_i_letters=0,
    fontsize=20,
    structure_size = 280,
    plots_size = 250,
    nflts = 1,
    ylims::Tuple = (:auto,nothing,nothing) # use (:manual,v1,v2) alternatively
  )
  @assert all(data.Objective .== "dh") "Expected objective to dh"
  IDs = data.ID
  ID_to_fig = Dict{String,Figure}()
  letters = collect('a':'z')[1+skip_i_letters:4+skip_i_letters]
  _bold(x) = rich(x)#,font=:bold)
  for ID in IDs
    loc_data = data[data.ID .== ID,:]
    @assert size(loc_data,1) == 1 "We require unique idenitfiers (IDs) for each optimisation solution"
    img_SOLID = load(joinpath(structure_path,first(loc_data.ID)*"_STRUCTURE.png"));
    img_VOID = load(joinpath(structure_path,first(loc_data.ID)*"_STRUCTURE_VOID.png"));
    with_theme(theme_latexfonts(), palette=(color=colors,markercolor=color,patchcolor=color)) do
      f = Figure(fontsize=fontsize, markersize = 24)
      if show_initial_struc
        axes = Axis(f[1,0],title = "Initial stucture",titlecolor=(:black,0.0),aspect = DataAspect(),yticksvisible=false,xticksvisible=false, width = structure_size, height = structure_size)
        Label(f[1, 0, TopLeft()], _bold("($(letters[1]))"),color=(:black,show_labels ? 1.0 : 0.0));
        img_init = load(joinpath(structure_path,loc_data.Initial_struc[1]*"_STRUCTURE.png"));
        image!(axes, rotr90(img_init))
        hidedecorations!(axes)
        hidespines!(axes)
      end
      axes = Axis(f[1,1],title = "Final Solid",titlecolor=(:black,0.0),aspect = DataAspect(),yticksvisible=false,xticksvisible=false, width = structure_size, height = structure_size)
      Label(f[1, 1, TopLeft()], (show_initial_struc ? _bold("($(letters[2]))") : _bold("($(letters[1]))")),color=(:black, show_labels && (~left_label_only || ~show_initial_struc) ? 1.0 : 0.0));
      image!(axes, rotr90(img_SOLID))
      hidedecorations!(axes)
      hidespines!(axes)
      axes = Axis(f[1,2],title = "Final void",titlecolor=(:black,0.0),aspect = DataAspect(),yticksvisible=false,xticksvisible=false, width = structure_size, height = structure_size)
      Label(f[1, 2, TopLeft()], (show_initial_struc ? _bold("($(letters[3]))") : _bold("($(letters[2]))")),color=(:black,show_labels && ~left_label_only ? 1.0 : 0.0));
      image!(axes, rotr90(img_VOID))
      hidedecorations!(axes)
      hidespines!(axes)
      req_c33 = cfmt("%.$(nflts)f",loc_data.Req_C33[1]*10); # Change to 10^9 units
      ax1 = Axis(f[1, 3], yticklabelcolor = colors[1], ylabelcolor = colors[1], xlabel = L"\text{Iteration}", ylabel=L"$\bar{d}_h(\Omega)\times10^{9}$",
        xtickalign=1,ytickalign=1,xminorticksvisible=true,yminorticksvisible=true,xgridvisible=true,ygridvisible=true,
        xminortickalign=1,yminortickalign=1,xminorticks=IntervalsBetween(5),yminorticks=IntervalsBetween(5),ytickformat = "{:.2f}", width = plots_size, height = plots_size)
      if ylims[1] == :manual
        ylims!(ax1,ylims[2],ylims[3])
      end
      ax2 = Axis(f[1, 3], yticklabelcolor = colors[2], ylabelcolor = colors[2], yaxisposition = :right, ylabel=L"""10^{-9}\bar{C}_{zzzz}(\Omega)-%$(req_c33)""",
        xtickalign=1,ytickalign=1,xminorticksvisible=true,yminorticksvisible=true,xgridvisible=true,ygridvisible=false,
        xminortickalign=1,yminortickalign=1,xminorticks=IntervalsBetween(5),yminorticks=IntervalsBetween(5),ytickformat = "{:.2f}", width = plots_size, height = plots_size,)
      hidespines!(ax2)
      hidexdecorations!(ax2)
      lines!(ax1,-first(loc_data.History).J,color=colors[1],linewidth = 3)
      lines!(ax2,first(loc_data.History).C33,linestyle=:dot,color=colors[2],linewidth = 3,label="constraint")
      Label(f[1, 2, TopRight()], show_initial_struc ? _bold("($(letters[4]))") : _bold("($(letters[3]))"),color=(:black,show_labels && ~left_label_only ? 1.0 : 0.0))#,padding=(-100,0,0,0))
      # Legend(f[0,3],ax2,orientation = :horizontal,margin=(0,0,-70,0))
      resize_to_layout!(f)
      ID_to_fig[ID]=f
    end
  end
  data[:,new_column_name] = getindex.((ID_to_fig,),data.ID);
  return data
end

function build_fig7(
  ordered_data;
  structure_path="$(@__DIR__)/../data/structures/raw/",
  colors=ColorSchemes.seaborn_deep.colors,
  show_labels = true,
  left_label_only = false,
  skip_i_letters=0,
  fontsize=20,
  structure_size = 280,
  plots_size = 250,
  nflts = 1)

  @assert length(ordered_data) == 3
  letters = collect('a':'z')[1+skip_i_letters:4+skip_i_letters]
  _bold(x) = rich(x)#,font=:bold)
  img_SOLID = map(d->load(joinpath(structure_path,d.ID*"_STRUCTURE.png")),ordered_data);
  with_theme(theme_latexfonts(), palette=(color=colors,markercolor=color,patchcolor=color)) do
    f = Figure(fontsize=fontsize, markersize = 24)
    function gen_fig!(j,i,frame_i,hist)
      axes = Axis(f[j,frame_i],titlecolor=(:black,0.0),aspect = DataAspect(),yticksvisible=false,xticksvisible=false, width = structure_size, height = structure_size)
      Label(f[j, frame_i, TopLeft()], _bold("($(letters[i]))"),color=(:black, show_labels && (~left_label_only || ~show_initial_struc) ? 1.0 : 0.0));
      image!(axes, rotr90(img_SOLID[i]))
      hidedecorations!(axes)
      hidespines!(axes)
      if hist
        # hist
        loc_data = ordered_data[i]
        req_c33 = cfmt("%.$(nflts)f",loc_data.Req_C33[1]*10); # Change to 10^9 units
        ax1 = Axis(f[j, frame_i+1], yticklabelcolor = colors[1], ylabelcolor = colors[1], xlabel = L"\text{Iteration}", ylabel=L"$\bar{d}_h(\Omega)\times10^{9}$",
          xtickalign=1,ytickalign=1,xminorticksvisible=true,yminorticksvisible=true,xgridvisible=true,ygridvisible=true,
          xminortickalign=1,yminortickalign=1,xminorticks=IntervalsBetween(5),yminorticks=IntervalsBetween(5),ytickformat = "{:.2f}", width = plots_size, height = plots_size)
        ax2 = Axis(f[j, frame_i+1], yticklabelcolor = colors[2], ylabelcolor = colors[2], yaxisposition = :right, ylabel=L"""10^{-9}\bar{C}_{zzzz}(\Omega)-%$(req_c33)""",
          xtickalign=1,ytickalign=1,xminorticksvisible=true,yminorticksvisible=true,xgridvisible=true,ygridvisible=false,
          xminortickalign=1,yminortickalign=1,xminorticks=IntervalsBetween(5),yminorticks=IntervalsBetween(5),ytickformat = "{:.2f}", width = plots_size, height = plots_size,)
        hidespines!(ax2)
        hidexdecorations!(ax2)
        lines!(ax1,-loc_data.History.J,color=colors[1],linewidth = 3)
        lines!(ax2,loc_data.History.C33,linestyle=:dot,color=colors[2],linewidth = 3,label="constraint")
      end
    end
    gen_fig!(1,1,1,true)
    gen_fig!(1,2,3,false)
    gen_fig!(1,3,4,false)
    resize_to_layout!(f)
    f
  end
end