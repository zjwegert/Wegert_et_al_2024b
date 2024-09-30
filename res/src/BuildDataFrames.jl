function build_max_bh_dataframe(;path=["/media/zachary/DATA/Data/Piezo/BVoigt/","/media/zachary/DATA/Data/Piezo/BReuss/"])
  dirs = readdir.(path)
  path_to_dir = Dict([path[i]=>dirs[i] for i = 1:length(path)]...);

  # Get data
  IDENTIFIER = String[];
  OBJ_NAME = String[];
  REQ_VOL_VAL = Float64[];
  REQ_DH_VAL = Float64[];
  FINAL_BH_VOIGT_VAL = Float64[];
  FINAL_BH_REUSS_VAL = Float64[];
  FINAL_VOL_VAL = Float64[];
  FINAL_DH_VAL = Float64[];
  FINAL_OBJ_CHNG = Vector{Float64}[];
  INIT_STRUC = String[];
  HISTORY_DATA = DataFrame[];
  EL = Int[];
  ORDER = Int[];
  GAMMA = Float64[];
  ALPHA_MIN = Float64[];
  AL_EXT_COEF = Float64[];
  LINE_SEARCH = Bool[];
  STIFFNESS_TENSOR = Matrix[];
  PIEZO_TENSOR = Matrix[];
  DIELECTRIC_TENSOR = Matrix[];

  for root_dir in keys(path_to_dir)
    for folder in path_to_dir[root_dir]
      if isdir(root_dir*folder)
        reg = match(r"3d_n=([0-9]+)_o=(\d)_max_([A-Za-z]+)_st_Vol=([0-9]*\.[0-9]+)_dh=([0-9]*\.[0-9]+)_gam=([0-9]*\.[0-9]+)_amin=([0-9]*\.[0-9]+)_alph=([0-9]*\.[0-9]+)_ls=([A-Za-z]+)_Struc=(.*)",folder)
        el = parse(Int,reg[1]); order = parse(Int,reg[2]); obj_name = reg[3]; vol_required = parse(Float64,reg[4]); dh_required = parse(Float64,reg[5]);
        gamma = parse(Float64,reg[6]); amin = parse(Float64,reg[7]); alpha_ext_coeff = parse(Float64,reg[8]); line_search = parse(Bool,reg[9]); init_struc = reg[10]
        history_data, history_header = readdlm(root_dir*folder*"/history.txt",',';header=true)
        history_dataframe = DataFrame(strip(history_header[1])=>history_data[:,1],strip(history_header[2])=>history_data[:,2],
          strip(history_header[3])=>history_data[:,3],strip(history_header[4])=>history_data[:,4])
        material_data = JSON3.read(root_dir*folder*"/material_data.json")
        Cᴴᵣₛ = reshape(material_data.Ch,6,6)
        Sᴴ = inv(Cᴴᵣₛ);
        BHₕ_Voigt = 1/9*(Cᴴᵣₛ[1,1]+Cᴴᵣₛ[2,2]+Cᴴᵣₛ[3,3]+2*(Cᴴᵣₛ[1,2]+Cᴴᵣₛ[1,3]+Cᴴᵣₛ[2,3]))
        Bᴴₕ_Reuss = (1/sum(Sᴴ[i,j] for i = 1:3, j = 1:3))

        _obj_name = obj_name=="Bh" ? "Bh_Reuss" : "Bh_Voigt"

        push!(IDENTIFIER,folder)
        push!(OBJ_NAME,_obj_name)
        push!(REQ_VOL_VAL,vol_required)
        push!(REQ_DH_VAL,dh_required)
        push!(FINAL_BH_VOIGT_VAL,BHₕ_Voigt)
        push!(FINAL_BH_REUSS_VAL,Bᴴₕ_Reuss)
        push!(FINAL_VOL_VAL,history_dataframe.Vol[end] + vol_required)
        push!(FINAL_DH_VAL,history_dataframe.dh[end] + dh_required)
        push!(FINAL_OBJ_CHNG, @.(abs(history_dataframe.J[end] - history_dataframe.J[end-5:end])/abs(history_dataframe.J[end])))
        push!(INIT_STRUC,init_struc)
        push!(HISTORY_DATA,history_dataframe)
        push!(EL,el)
        push!(ORDER,order)
        push!(GAMMA,gamma)
        push!(ALPHA_MIN,amin)
        push!(AL_EXT_COEF,alpha_ext_coeff)
        push!(LINE_SEARCH,line_search)
        push!(STIFFNESS_TENSOR,reshape(material_data.Ch,6,6))
        push!(PIEZO_TENSOR,reshape(material_data.eh,3,6))
        push!(DIELECTRIC_TENSOR,reshape(material_data.Kh,3,3))
      end
    end
  end

  data=DataFrame(ID=IDENTIFIER,Objective=OBJ_NAME,Req_Vol=REQ_VOL_VAL,Req_dh=REQ_DH_VAL,
    Bh_Voigt=FINAL_BH_VOIGT_VAL,Bh_Reuss=FINAL_BH_REUSS_VAL,Vol=FINAL_VOL_VAL,
    dh=FINAL_DH_VAL,ObjChange=FINAL_OBJ_CHNG,Initial_struc=INIT_STRUC,History=HISTORY_DATA,
    C=STIFFNESS_TENSOR,e=PIEZO_TENSOR,K=DIELECTRIC_TENSOR,
    n=EL,o=ORDER,γ=GAMMA,α_min=ALPHA_MIN,α=AL_EXT_COEF,line_search=LINE_SEARCH)

  jldsave("$(@__DIR__)/../data/max_bh_results.jld2";data)
end

function build_max_dh_dataframe(;path="/media/zachary/DATA/Data/Piezo/dh/")
  dirs = readdir.(path)
  path_to_dir = Dict([path[i]=>dirs[i] for i = 1:length(path)]...);

  # Get data
  IDENTIFIER = String[];
  OBJ_NAME = String[];
  REQ_C33_VAL = Float64[];
  FINAL_BH_VOIGT_VAL = Float64[];
  FINAL_BH_REUSS_VAL = Float64[];
  FINAL_VOL_VAL = Float64[];
  FINAL_DH_VAL = Float64[];
  FINAL_OBJ_CHNG = Vector{Float64}[];
  INIT_STRUC = String[];
  HISTORY_DATA = DataFrame[];
  EL = Int[];
  ORDER = Int[];
  GAMMA = Float64[];
  ALPHA_MIN = Float64[];
  AL_EXT_COEF = Float64[];
  LINE_SEARCH = Bool[];
  STIFFNESS_TENSOR = Matrix[];
  PIEZO_TENSOR = Matrix[];
  DIELECTRIC_TENSOR = Matrix[];

  for root_dir in keys(path_to_dir)
    for folder in path_to_dir[root_dir]
      if isdir(root_dir*folder)
        reg = match(r"3d_n=([0-9]+)_o=(\d)_max_([A-Za-z]+)_st_C33=([0-9]*\.[0-9]+)_gam=([0-9]*\.[0-9]+)_amin=([0-9]*\.[0-9]+)_alph=([0-9]*\.[0-9]+)_ls=([A-Za-z]+)_Struc=(.*)",folder)
        el = parse(Int,reg[1]); order = parse(Int,reg[2]); obj_name = reg[3]; c33_required = parse(Float64,reg[4]);
        gamma = parse(Float64,reg[5]); amin = parse(Float64,reg[6]); alpha_ext_coeff = parse(Float64,reg[7]); line_search = parse(Bool,reg[8]); init_struc = reg[9]
        history_data, history_header = readdlm(root_dir*folder*"/history.txt",',';header=true)
        history_dataframe = DataFrame(strip(history_header[1])=>history_data[:,1],strip(history_header[2])=>history_data[:,2],
          strip(history_header[3])=>history_data[:,3])
        material_data = JSON3.read(root_dir*folder*"/material_data.json")
        Cᴴᵣₛ = reshape(material_data.Ch,6,6)
        Sᴴ = inv(Cᴴᵣₛ);
        BHₕ_Voigt = 1/9*(Cᴴᵣₛ[1,1]+Cᴴᵣₛ[2,2]+Cᴴᵣₛ[3,3]+2*(Cᴴᵣₛ[1,2]+Cᴴᵣₛ[1,3]+Cᴴᵣₛ[2,3]))
        Bᴴₕ_Reuss = (1/sum(Sᴴ[i,j] for i = 1:3, j = 1:3))

        push!(IDENTIFIER,folder)
        push!(OBJ_NAME,obj_name)
        push!(REQ_C33_VAL,c33_required)
        push!(FINAL_BH_VOIGT_VAL,BHₕ_Voigt)
        push!(FINAL_BH_REUSS_VAL,Bᴴₕ_Reuss)
        push!(FINAL_VOL_VAL,material_data.Vol)
        push!(FINAL_DH_VAL,material_data.dh)
        push!(FINAL_OBJ_CHNG, @.(abs(history_dataframe.J[end] - history_dataframe.J[end-5:end])/abs(history_dataframe.J[end])))
        push!(INIT_STRUC,init_struc)
        push!(HISTORY_DATA,history_dataframe)
        push!(EL,el)
        push!(ORDER,order)
        push!(GAMMA,gamma)
        push!(ALPHA_MIN,amin)
        push!(AL_EXT_COEF,alpha_ext_coeff)
        push!(LINE_SEARCH,line_search)
        push!(STIFFNESS_TENSOR,reshape(material_data.Ch,6,6))
        push!(PIEZO_TENSOR,reshape(material_data.eh,3,6))
        push!(DIELECTRIC_TENSOR,reshape(material_data.Kh,3,3))
      end
    end
  end

  data=DataFrame(ID=IDENTIFIER,Objective=OBJ_NAME,Req_C33=REQ_C33_VAL,
    Bh_Voigt=FINAL_BH_VOIGT_VAL,Bh_Reuss=FINAL_BH_REUSS_VAL,Vol=FINAL_VOL_VAL,
    dh=FINAL_DH_VAL,ObjChange=FINAL_OBJ_CHNG,Initial_struc=INIT_STRUC,History=HISTORY_DATA,
    C=STIFFNESS_TENSOR,e=PIEZO_TENSOR,K=DIELECTRIC_TENSOR,
    n=EL,o=ORDER,γ=GAMMA,α_min=ALPHA_MIN,α=AL_EXT_COEF,line_search=LINE_SEARCH)

  jldsave("$(@__DIR__)/../data/max_dh_results.jld2";data)
end

function build_benchmark_dataframe(;path="$(@__DIR__)/../data/solver_benchmarks/")
  files = readdir(path)

  NCPU = Int[];
  N = Int[];
  SOLVER = String[];
  RTOL = Float64[];
  TIME = Float64[];
  U_ITS = Int[];
  PHI_ITS = Int[];
  OUTER_ITS = Int[];
  RES = Float64[];

  for file in files
    if isfile(path*file)
      reg = match(r"3d_ncpus=([0-9]+)_n=([0-9]+)_o=1_solver=([A-Za-z]+)_Urtol=(1e-[0-9]+)_Phirtol=(1e-[0-9]+).txt",file)
      ncpu = reg[1]; n = reg[2]; solver = reg[3]; tol = reg[4];
      data = readdlm(path*file,',',skipstart=1)

      push!(NCPU,parse(Int,ncpu)); push!(N,parse(Int,n)); push!(SOLVER,solver); push!(RTOL,parse(Float64,tol));
      push!(TIME,maximum(data[:,1])); push!(U_ITS,convert(Int64,maximum(data[:,2])));
      push!(PHI_ITS,convert(Int64,maximum(data[:,3]))); push!(OUTER_ITS,convert(Int64,maximum(data[:,4])));
      push!(RES,maximum(data[:,5]));
    end
  end

  data = DataFrame(N=NCPU,n=N,Solver=SOLVER,rtol=RTOL,time=TIME,U_its=U_ITS,phi_its=PHI_ITS,outer_its=OUTER_ITS,residual=RES)
  jldsave("$(@__DIR__)/../data/benchmark_data.jld2";data)
end