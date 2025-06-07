### Post-processing of microstructure results using CutFEM.

## Run on most recent versions of Gridap, GridapTopOpt, GridapEmbedded, GridapDistributed, etc.
# Notes:
#  - Need to use a dev version of GridapTopOpt and disable the compats (Project.toml)
#  - Need to disable compats for this package (Project.toml)

using Wegert_et_al_2024b,GridapTopOpt,Gridap,Gridap.TensorValues,JSON3,Comonicon
# Distributed Dependencies
using GridapDistributed, GridapPETSc, GridapSolvers, PartitionedArrays,
  SparseMatricesCSR, GridapSolvers.BlockSolvers, Gridap.MultiField

using Gridap.FESpaces
using GridapEmbedded, GridapEmbedded.LevelSetCutters
using GridapDistributed: allocate_in_range, allocate_in_domain

function run(ranks,params)
  mesh_partition,n,path,struc_name = params

  order = 1
  el_size = (n,n,n)
  xmax,ymax,zmax=(1.0,1.0,1.0)
  dom = (0,xmax,0,ymax,0,zmax)

  ## FE Setup
  model = CartesianDiscreteModel(ranks,mesh_partition,dom,el_size,isperiodic=(true,true,true))
  f_Γ_D(x) = iszero(x)
  update_labels!(1,model,f_Γ_D,"origin")

  ## Level-set function
  reffe_scalar = ReferenceFE(lagrangian,Float64,order)
  V_φ = TestFESpace(model,reffe_scalar)
  φh = interpolate(0,V_φ)
  pload!(joinpath(path,struc_name),get_free_dof_values(φh))
  consistent!(get_free_dof_values(φh)) |> fetch

  ## Triangulations and measures
  Ω = Triangulation(model)
  dΩ = Measure(Ω,2*order)
  vol_D = sum(∫(1)dΩ)

  geo = DiscreteGeometry(φh,model)
  cutgeo = cut(model,geo)
  Ω_act = Triangulation(cutgeo,ACTIVE)
  Ω_phys = Triangulation(cutgeo,PHYSICAL)
  dΩ_phys = Measure(Ω_phys,2*order)

  Γg = GhostSkeleton(cutgeo)
  dΓg = Measure(Γg,2*order)
  n_Γg = get_normal_vector(Γg)

  ## Spaces
  reffe = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
  V = TestFESpace(Ω_act,reffe;conformity=:H1,dirichlet_tags=["origin"])
  U = TrialFESpace(V,VectorValue(0.0,0.0,0.0))
  Q = TestFESpace(Ω_act,reffe_scalar;conformity=:H1,dirichlet_tags=["origin"])
  P = TrialFESpace(Q,0)
  UP = MultiFieldFESpace([U,P])
  VQ = MultiFieldFESpace([V,Q])

  ## Operators
  mat = PZT5A(3)
  C, e, κ, A₁₁, A₁₂, A₂₂, _, _, _ = mat
  εᴹ = (SymTensorValue(1.0,0.0,0.0,0.0,0.0,0.0),
        SymTensorValue(0.0,0.0,0.0,1.0,0.0,0.0),
        SymTensorValue(0.0,0.0,0.0,0.0,0.0,1.0),
        SymTensorValue(0.0,0.0,0.0,0.0,0.5,0.0),
        SymTensorValue(0.0,0.0,0.5,0.0,0.0,0.0),
        SymTensorValue(0.0,0.5,0.0,0.0,0.0,0.0))
  Eⁱ = (VectorValue(1.,0.,0.),VectorValue(0.,1.,0),VectorValue(0.,0.,1.))

  h = 1/100
  ## Ghost penalty for the potential, based on Section 3.2 of 10.1016/j.cma.2017.09.005
  # 10 is a scaling factor based on the rough size of largest commonent in normalised stiffness tensor
  α_Gd = 0.1
  γ_Gd = α_Gd*10*h^3
  j_u(d,s) = γ_Gd*(jump(n_Γg ⋅ ∇(s)) ⋅ jump(n_Γg ⋅ ∇(d)))
  ## Ghost penalty for the potential, based on Section 6.1 of 10.1002/nme.4823
  # 40 is a scaling factor based on the rough size of largest commonent in normalised dielectric tensor
  γ_Gd_2 = γ_Gd*40*h
  j_ϕ(ϕ,q) = γ_Gd_2*jump(n_Γg ⋅ ∇(q))*jump(n_Γg ⋅ ∇(ϕ))

  a((u,ϕ),(v,q)) = ∫((A₁₁*((C ⊙ ε(u)) ⊙ ε(v)) - A₁₂*((-∇(ϕ) ⋅ e) ⊙ ε(v)) +
                      A₁₂*((e ⋅² ε(u)) ⋅ -∇(q)) + A₂₂*((κ ⋅ -∇(ϕ)) ⋅ -∇(q))) )dΩ_phys +
                    ∫(j_u(u,v) + j_ϕ(ϕ,q))dΓg;
  l_ε = [((v,q),) -> ∫((-A₁₁*((C ⊙ εᴹ[i]) ⊙ ε(v)) - A₁₂*((e ⋅² εᴹ[i]) ⋅ -∇(q))) )dΩ_phys for i = 1:length(εᴹ)]
  l_E = [((v,q),) -> ∫(( A₁₂*((Eⁱ[i] ⋅ e) ⊙ ε(v))  - A₂₂*((κ ⋅ Eⁱ[i]) ⋅ -∇(q))) )dΩ_phys for i = 1:length(Eⁱ)]
  li = [l_ε; l_E]

  ## Solver
  solver = PETScLinearSolver()
  assem = SparseMatrixAssembler(UP,VQ)
  K = assemble_matrix(a,assem,UP,VQ)
  ns = numerical_setup(symbolic_setup(solver,K),K)
  numerical_setup!(ns,K)

  bi = [allocate_in_range(K) for _ = 1:length(li)];
  xi = [allocate_in_domain(K) for _ = 1:length(li)];

  for i in eachindex(li)
    assemble_vector!(li[i],bi[i],assem,VQ)
    solve!(xi[i],ns,bi[i])
  end

  xhi = [FEFunction(UP,xi[i]) for i = 1:length(li)];
  uφ = collect(Iterators.flatten(xhi))

  ## Functionals
  N = length(εᴹ)
  function Cᴴ(r,s,uϕ)
      u_s = uϕ[2s-1]; ϕ_s = uϕ[2s]
      ∫(A₁₁*((C ⊙ (ε(u_s) + εᴹ[s])) ⊙ εᴹ[r]) - A₁₂*((-∇(ϕ_s) ⋅ e) ⊙ εᴹ[r]))dΩ_phys;
  end
  function eᴴ(i,s,uϕ)
    u_s = uϕ[2s-1]; ϕ_s = uϕ[2s]
    ∫(A₁₂*((e ⋅² (ε(u_s) + εᴹ[s])) ⋅ Eⁱ[i]) + A₂₂*((κ ⋅ (-∇(ϕ_s))) ⋅ Eⁱ[i]))dΩ_phys;
  end
  function κᴴ(i,j,uϕ)
    u_j = uϕ[2N+2j-1]; ϕ_j = uϕ[2N+2j]
    ∫(A₁₂*((e ⋅² ε(u_j)) ⋅ Eⁱ[i]) + A₂₂*((κ ⋅ (-∇(ϕ_j) + Eⁱ[j])) ⋅ Eⁱ[i]))dΩ_phys;
  end

  ## Write material values
  Cᴴᵣₛ = map(rs->sum(Cᴴ(rs[1],rs[2],uφ)), CartesianIndices((1:6, 1:6)))*mat.c1/mat.A₁₁
  eᴴᵢₛ = map(rs->sum(eᴴ(rs[1],rs[2],uφ)), CartesianIndices((1:3, 1:6)))*mat.e1/mat.A₁₂
  κᴴᵢⱼ = map(rs->sum(κᴴ(rs[1],rs[2],uφ)), CartesianIndices((1:3, 1:3)))*mat.κ1/mat.A₂₂
  Sᴴ = inv(Cᴴᵣₛ)
  dᴴ = eᴴᵢₛ*Sᴴ;
  dᴴₕ = (dᴴ[3,1] + dᴴ[3,2] + dᴴ[3,3])*(mat.e1/mat.A₁₂)
  Bᴴₕ_Voigt = 1/9*(Cᴴᵣₛ[1,1]+Cᴴᵣₛ[2,2]+Cᴴᵣₛ[3,3]+2*(Cᴴᵣₛ[1,2]+Cᴴᵣₛ[1,3]+Cᴴᵣₛ[2,3]))
  Bᴴₕ_Reuss = (1/sum(Sᴴ[i,j] for i = 1:3, j = 1:3))
  Vol = sum(∫(1/vol_D)dΩ_phys);
  data = Dict("Ch"=>Cᴴᵣₛ,"eh"=>eᴴᵢₛ,"Kh"=>κᴴᵢⱼ,"dh"=>dᴴₕ,"Bh_Voigt"=>Bᴴₕ_Voigt,"Bh_Reuss"=>Bᴴₕ_Reuss,"Vol"=>Vol)
  i_am_main(ranks) && open(joinpath(path,"material_data_aggfem.json"),"w") do io
    JSON3.pretty(io, data)
  end
end

function main(px::Int,py::Int,pz::Int,n::Int,path::String)
  mesh_partition = (px,py,pz)
  with_mpi() do distribute
    ranks = distribute(LinearIndices((prod(mesh_partition),)))
    struc_name = last(readdir(path))
    if i_am_main(ranks)
      println("--------------------------------------------------")
      println("Parsed args for post-processing")
      println("       (px,py,pz) => ", mesh_partition)
      println("                n => ", n)
      println("             path => ", path)
      println("       struc_name => ", struc_name)
      println("--------------------------------------------------")
    end
    lsf_fn(x) = invokelatest(eval(Meta.parse(lsf_func)),x)
    params = (;mesh_partition,n,path,struc_name)
    petsc_options = "-ksp_converged_reason -ksp_error_if_not_converged true -pc_type lu -pc_factor_mat_solver_type superlu_dist"
    GridapPETSc.with(;args=split(petsc_options)) do
      run(ranks,params)
    end
  end
end

main(parse(Int,ARGS[1]),parse(Int,ARGS[2]),parse(Int,ARGS[3]),parse(Int,ARGS[4]),ARGS[5])