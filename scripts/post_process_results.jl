using Wegert_et_al_2024b,GridapTopOpt,Gridap,Gridap.TensorValues,GridapEmbedded,JSON3,Comonicon
# Distributed Dependencies
using GridapDistributed, GridapPETSc, GridapSolvers, PartitionedArrays,
  SparseMatricesCSR, GridapSolvers.BlockSolvers, Gridap.MultiField

function main(ranks,params)
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

  ## Spaces
  strategy = AggregateAllCutCells()
  aggregates = aggregate(strategy,cutgeo)

  reffe = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
  Vstd = TestFESpace(Ω_act,reffe;conformity=:H1,dirichlet_tags=["origin"])
  V = AgFEMSpace(Vstd,aggregates)
  U = TrialFESpace(V,VectorValue(0.0,0.0,0.0))
  Qstd = TestFESpace(Ω_act,reffe_scalar;conformity=:H1,dirichlet_tags=["origin"])
  Q = AgFEMSpace(Qstd,aggregates)
  P = TrialFESpace(Q,0)
  mfs = BlockMultiFieldStyle()
  UP = MultiFieldFESpace([U,P];style=mfs)
  VQ = MultiFieldFESpace([V,Q];style=mfs)

  ## Operators
  C, e, κ, A₁₁, A₁₂, A₂₂, _, _, _ = PZT5A(3)
  εᴹ = (SymTensorValue(1.0,0.0,0.0,0.0,0.0,0.0),
        SymTensorValue(0.0,0.0,0.0,1.0,0.0,0.0),
        SymTensorValue(0.0,0.0,0.0,0.0,0.0,1.0),
        SymTensorValue(0.0,0.0,0.0,0.0,0.5,0.0),
        SymTensorValue(0.0,0.0,0.5,0.0,0.0,0.0),
        SymTensorValue(0.0,0.5,0.0,0.0,0.0,0.0))
  Eⁱ = (VectorValue(1.,0.,0.),VectorValue(0.,1.,0),VectorValue(0.,0.,1.))
  a((u,ϕ),(v,q)) = ∫((A₁₁*((C ⊙ ε(u)) ⊙ ε(v)) - A₁₂*((-∇(ϕ) ⋅ e) ⊙ ε(v)) +
                      A₁₂*((e ⋅² ε(u)) ⋅ -∇(q))  + A₂₂*((κ ⋅ -∇(ϕ)) ⋅ -∇(q))) )dΩ_phys;
  l_ε = [((v,q)) -> ∫((-A₁₁*((C ⊙ εᴹ[i]) ⊙ ε(v)) - A₁₂*((e ⋅² εᴹ[i]) ⋅ -∇(q))) )dΩ_phys for i = 1:length(εᴹ)]
  l_E = [((v,q)) -> ∫(( A₁₂*((Eⁱ[i] ⋅ e) ⊙ ε(v))  - A₂₂*((κ ⋅ Eⁱ[i]) ⋅ -∇(q))) )dΩ_phys for i = 1:length(Eⁱ)]
  li = [l_ε; l_E]

  ## Solver
  Tm = SparseMatrixCSR{0,PetscScalar,PetscInt}
  Tv = Vector{PetscScalar}
  solver_u = ElasticitySolver(V;rtol=10^-2)
  solver_ϕ = PETScLinearSolver(cgamg_ksp_setup(rtol=10^-2))

  P_solver = GridapSolvers.BlockTriangularSolver([solver_u,solver_ϕ],half=:lower);
  schur_preconditioner = BlockSchurPreconditioner(P_solver,a,UP,VQ,Ω_act,10^-10;
    schur_assembler = SparseMatrixAssembler(Tm,Tv,P,Q))
  solver = GridapSolvers.LinearSolvers.GMRESSolver(100;Pr=schur_preconditioner,rtol=1.e-12,
    verbose=i_am_main(ranks))

  assem = SparseMatrixAssembler(UP,VQ)
  K = assemble_matrix(a_fwd,assem,UP,VQ)
  ns = numerical_setup(symbolic_setup(solver,K),K)
  numerical_setup!(ns,K)

  bi = [allocate_in_range(K) for _ = 1:length(li)];
  xi = [allocate_in_domain(K) for _ = 1:length(li)];

  for i in eachindex(li)
    assemblze_vector!(li[i],bi[i],assem,VQ)
    solve!(xi[i],ns,bi[i])
  end

  xhi = [FEFunction(UP,xi[i]) for i = 1:length(li)];
  uφ = collect(Iterators.flatten(xhi))

  ## Functionals
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
  dᴴₕ = (dᴴ[3,1] + dᴴ[3,2] + dᴴ[3,3])*(mat.e1/mat.A₁₂)*(mat.c1/mat.A₁₁)^-1
  Bᴴₕ_Voigt = 1/9*(Cᴴᵣₛ[1,1]+Cᴴᵣₛ[2,2]+Cᴴᵣₛ[3,3]+2*(Cᴴᵣₛ[1,2]+Cᴴᵣₛ[1,3]+Cᴴᵣₛ[2,3]))
  Bᴴₕ_Reuss = (1/sum(Sᴴ[i,j] for i = 1:3, j = 1:3))
  Vol = sum(∫(1/vol_D)dΩ_phys);
  data = Dict("Ch"=>Cᴴᵣₛ,"eh"=>eᴴᵢₛ,"Kh"=>κᴴᵢⱼ,"dh"=>dᴴₕ,"Bh_Voigt"=>Bᴴₕ_Voigt,"Bh_Reuss"=>Bᴴₕ_Reuss,"Vol"=>Vol)
  i_am_main(ranks) && open(joinpath(path,"material_data_aggfem.json"),"w") do io
    JSON3.pretty(io, data)
  end
end

"""
Run the above post-processing solve.

# Arguments

# Options

- `--px <arg>`: partition along x-axis.
- `--py <arg>`: partition along y-axis.
- `--pz <arg>`: partition along z-axis.
- `--n <arg>` : mesh partition along each axis.
- `--path <arg>`: path to the directory where the results are stored.

# Flags

"""
@main function run(;
    px::Int,
    py::Int,
    pz::Int,
    n::Int=100,
    path::String)
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
    GridapPETSc.with(args=split(petsc_opts)) do
      main(ranks,params)
    end
  end
end