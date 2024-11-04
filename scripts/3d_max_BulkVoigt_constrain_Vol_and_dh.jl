using Wegert_et_al_2024b,GridapTopOpt,Gridap,Gridap.TensorValues,JSON3,Comonicon
# Distributed Dependencies
using GridapDistributed, GridapPETSc, GridapSolvers, PartitionedArrays,
  SparseMatricesCSR, GridapSolvers.BlockSolvers, Gridap.MultiField

using GridapTopOpt_Piezo: bulk_modulus_voigt_3d, hydrostatic_coupling_3d

function main(ranks,params)
  mesh_partition,dh,vf,n,order,γ,amin,amax,alpha_coeff,lsf_fn,path,
    iter_mod,solver_method,petsc_opts,init_struc,reinit_mod,
    line_search,xi,xi_reduce_coef,xi_reduce_abs_tol,elast_solver_rtol,cond_solver_rtol = params

  el_size = (n,n,n)
  xmax,ymax,zmax=(1.0,1.0,1.0)
  dom = (0,xmax,0,ymax,0,zmax)
  γ_reinit = 0.5
  max_steps = floor(Int,order*minimum(el_size)/10)
  tol = 1/(5order^2)/minimum(el_size)
  η_coeff = 2
  α_coeff = alpha_coeff#*max_steps*γ
  jld2_path = path*"/lsf_files/"
  i_am_main(ranks) && mkpath(jld2_path)

  sym = isequal(solver_method,:Schur) ? false : true

  ## FE Setup
  model = CartesianDiscreteModel(ranks,mesh_partition,dom,el_size,isperiodic=(true,true,true))
  el_Δ = get_el_Δ(model)
  f_Γ_D(x) = iszero(x)
  update_labels!(1,model,f_Γ_D,"origin")

  ## Triangulations and measures
  Ω = Triangulation(model)
  dΩ = Measure(Ω,2*order)
  vol_D = sum(∫(1)dΩ)

  ## Spaces
  reffe = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
  reffe_scalar = ReferenceFE(lagrangian,Float64,order)
  V = TestFESpace(model,reffe;conformity=:H1,dirichlet_tags=["origin"])
  U = TrialFESpace(V,VectorValue(0.0,0.0,0.0))
  Q = TestFESpace(model,reffe_scalar;conformity=:H1,dirichlet_tags=["origin"])
  P = TrialFESpace(Q,0)
  mfs = BlockMultiFieldStyle()
  UP = MultiFieldFESpace([U,P];style=mfs)
  VQ = MultiFieldFESpace([V,Q];style=mfs)

  V_φ = TestFESpace(model,reffe_scalar)
  V_reg = TestFESpace(model,reffe_scalar)
  U_reg = TrialFESpace(V_reg)

  ## Create FE functions
  φh = interpolate(lsf_fn,V_φ)
  if ~isempty(init_struc)
    pload!(init_struc,get_free_dof_values(φh))
    consistent!(get_free_dof_values(φh)) |> fetch
  end

  ## Interpolation and weak form
  interp = SmoothErsatzMaterialInterpolation(η = η_coeff*maximum(el_Δ),ϵ=10^-6)
  H,DH,ρ = interp.H,interp.DH,interp.ρ

  ## Material tensors and Weak form
  mat = PZT5A(3)
  a,l,εᴹ,Eⁱ = piezo_hom_weak_form(mat,interp;sym)

  ## Optimisation functionals
  Cᴴ,eᴴ,κᴴ,DCᴴ,Deᴴ,Dκᴴ = homogenised_mat_tensor_functionals(mat,interp,εᴹ,Eⁱ)
  dᴴ,Ddᴴ = hydrostatic_coupling_3d(Cᴴ,eᴴ,κᴴ,DCᴴ,Deᴴ,Dκᴴ)
  Bᴴ,DBᴴ = bulk_modulus_voigt_3d(Cᴴ,eᴴ,κᴴ,DCᴴ,Deᴴ,Dκᴴ)

  # Problem
  J(uϕ,φ,dΩ) = -1*Bᴴ(uϕ,φ,dΩ)
  DJ(q,uϕ,φ,dΩ) = -1*DBᴴ(q,uϕ,φ,dΩ)
  C1(uϕ,φ,dΩ) = ∫(((ρ ∘ φ) - vf)/vol_D)dΩ;
  DC1(q,uϕ,φ,dΩ) = ∫(-1/vol_D*q*(DH ∘ φ)*(norm ∘ ∇(φ)))dΩ
  C2(uϕ,φ,dΩ) = dᴴ(uϕ,φ,dΩ) - ∫(dh/vol_D)dΩ
  DC2(q,uϕ,φ,dΩ) = Ddᴴ(q,uϕ,φ,dΩ)

  ## Finite difference solver and level set function
  stencil = HamiltonJacobiEvolution(FirstOrderStencil(3,Float64),model,V_φ,tol,max_steps)

  ## Setup solver and FE operators
  Tm = SparseMatrixCSR{0,PetscScalar,PetscInt}
  Tv = Vector{PetscScalar}
  solver_u = ElasticitySolver(V;rtol=elast_solver_rtol)
  solver_ϕ = PETScLinearSolver(cgamg_ksp_setup(rtol=cond_solver_rtol))

  if isequal(solver_method,:Schur)
    P_solver = GridapSolvers.BlockTriangularSolver([solver_u,solver_ϕ],half=:lower);
    schur_preconditioner = BlockSchurPreconditioner(P_solver,(u,v)->a(u,v,φh,dΩ),UP,VQ,Ω,10^-10;
      schur_assembler = SparseMatrixAssembler(Tm,Tv,P,Q))
    solver = GridapSolvers.LinearSolvers.GMRESSolver(100;Pr=schur_preconditioner,rtol=1.e-12,
      verbose=i_am_main(ranks))
  elseif isequal(solver_method,:TriCG)
    solver = TriCGSolver(solver_u,solver_ϕ;rtol=1.e-8,verbose=i_am_main(ranks))
  elseif isequal(solver_method,:BlockGMRES)
    P = GridapSolvers.BlockDiagonalSolver([solver_u,solver_ϕ])
    solver = GridapSolvers.LinearSolvers.GMRESSolver(100;Pr=P,rtol=1.e-8,verbose=i_am_main(ranks))
  else
    error("Solver ($solver_method) not defined")
  end

  state_map = RepeatingAffineFEStateMap(9,a,l,UP,VQ,V_φ,U_reg,φh,dΩ;
    assem_U = SparseMatrixAssembler(Tm,Tv,UP,VQ),
    assem_adjoint = SparseMatrixAssembler(Tm,Tv,VQ,UP),
    assem_deriv = SparseMatrixAssembler(Tm,Tv,U_reg,U_reg),
    ls = solver,adjoint_ls = solver
  )
  pcfs = PDEConstrainedFunctionals(J,[C1,C2],state_map;analytic_dJ=DJ,analytic_dC=[DC1,DC2])

  ## Hilbertian extension-regularisation problems
  α = α_coeff*maximum(el_Δ)
  a_hilb(p,q) = ∫(α^2*∇(p)⋅∇(q) + p*q)dΩ;
  vel_ext = VelocityExtension(
    a_hilb,U_reg,V_reg;
    assem = SparseMatrixAssembler(Tm,Tv,U_reg,V_reg),
    ls = PETScLinearSolver(cgamg_ksp_setup(rtol=10^-8))
  )

  ## Optimiser
  optimiser = HilbertianProjection(pcfs,stencil,vel_ext,φh;
    γ,γ_reinit,ls_γ_max=γ,α_min=amin,α_max=amax,ls_ξ=xi,ls_ξ_reduce_abs_tol=xi_reduce_abs_tol,
    ls_ξ_reduce_coef=xi_reduce_coef,ls_enabled=line_search,reinit_mod,
    verbose=i_am_main(ranks),constraint_names=[:Vol,:dh])
  for (it, uh, φh) in optimiser
    if iszero(it % iter_mod)
      psave(jld2_path*"out$it",get_free_dof_values(φh))
    end
    write_history(path*"/history.txt",optimiser.history)
  end
  it = get_history(optimiser).niter
  psave(jld2_path*"out$it",get_free_dof_values(φh))
  writevtk(Ω,path*"out$it",cellfields=["φ"=>φh,"H(φ)"=>(H ∘ φh),"|∇(φ)|"=>(norm ∘ ∇(φh))])

  # Write material values
  uh = get_state(state_map)
  Cᴴᵣₛ = map(rs->sum(Cᴴ(rs[1],rs[2],uh,φh,dΩ)), CartesianIndices((1:6, 1:6)))*mat.c1/mat.A₁₁
  eᴴᵢₛ = map(rs->sum(eᴴ(rs[1],rs[2],uh,φh,dΩ)), CartesianIndices((1:3, 1:6)))*mat.e1/mat.A₁₂
  κᴴᵢⱼ = map(rs->sum(κᴴ(rs[1],rs[2],uh,φh,dΩ)), CartesianIndices((1:3, 1:3)))*mat.κ1/mat.A₂₂
  Sᴴ = inv(Cᴴᵣₛ)
  dᴴₕ = sum(dᴴ(uh,φh,dΩ))*(mat.e1/mat.A₁₂)*(mat.c1/mat.A₁₁)^-1
  Bᴴₕ_Voigt = 1/9*(Cᴴᵣₛ[1,1]+Cᴴᵣₛ[2,2]+Cᴴᵣₛ[3,3]+2*(Cᴴᵣₛ[1,2]+Cᴴᵣₛ[1,3]+Cᴴᵣₛ[2,3]))
  Bᴴₕ_Reuss = (1/sum(Sᴴ[i,j] for i = 1:3, j = 1:3))
  Vol = sum(∫((ρ ∘ φh)/vol_D)dΩ);
  data = Dict("Ch"=>Cᴴᵣₛ,"eh"=>eᴴᵢₛ,"Kh"=>κᴴᵢⱼ,"dh"=>dᴴₕ,"Bh_Voigt"=>Bᴴₕ_Voigt,"Bh_Reuss"=>Bᴴₕ_Reuss,"Vol"=>Vol)
  i_am_main(ranks) && open(path*"material_data.json","w") do io
    JSON3.pretty(io, data)
  end
end

"""
Run the above optimisation problem. Note that the default value of
write_dir requires that the environment variable PROJECT_DIR is defined.

# Arguments

# Options

- `--px <arg>`: partition along x-axis.
- `--py <arg>`: partition along y-axis.
- `--pz <arg>`: partition along z-axis.
- `--n <arg>` : mesh partition along each axis.
- `--dh <arg>`: value of hydrostatic coupling constraint
- `--vol <arg>`: value of volume constraint
- `--amin <arg>`: value of alpha min in HPM
- `--amax <arg>`: value of alpha max in HPM
- `--xi <arg>`: value of ξ in line search test (J_new < J + ξ*abs(J))
- `--xi_reduce_coef <arg>`: Multiplier on ξ when constraints are within absolute tolerance defined by xi_reduce_abs_tol
- `--xi_reduce_abs_tol <arg>`: Tolerance on constraints to reduce ξ using xi_reduce_coef
- `--order <arg>`: FE order
- `--gamma <arg>`: coeffient on the time step size for solving the Hamilton-Jacobi evolution equation.
- `--alpha_coeff <arg>`: coeffient on alpha value
- `--lsf_func <arg>`: initial level set function. This gets parsed and evaluated so it should be something like
  initial_lsf(2,0.2) or an anonymous function x->-sqrt((x[1]-0.5)^2+(x[2]-0.5)^2)+0.55
- `--write_dir_suffix <arg>`: helper to add text to end of write dir (nothing by default)
- `--write_dir <arg>`: directory where results are written
- `--iter_mod <arg>`: how often results are written to write_dir
- `--reinit_mod <arg>`: how often to reinitialise
- `--solver_method <arg>`: solver to be used (Schur,TriCG,BlockGMRES)
- `--petsc_opts <arg>`: options for PETScLinearSolver
- `--init_struc <arg>`: directory for a starting structure, empty string by default
- `--elast_solver_rtol <arg>`: rtol for elasticity solver
- `--cond_solver_rtol <arg>`: rtol for conductivity solver

# Flags

- `--line_search`: Set if line search is used in HPM
"""
@main function run(;
    px::Int,
    py::Int,
    pz::Int,
    dh::Float64,
    vol::Float64,
    n::Int=100,
    order::Int=1,
    gamma::Float64=0.05,
    amin::Float64=0.5,
    amax::Float64=1.0,
    alpha_coeff::Float64=4.0,
    line_search::Bool=false,
    xi::Float64=1.0,
    xi_reduce_coef::Float64=0.0025,
    xi_reduce_abs_tol::Float64=0.01,
    lsf_func="initial_lsf(2,0.4)",
    write_dir_suffix::String="",
    write_dir::String=ENV["PROJECT_DIR"]*"/results/3d_n=$(n)_o=$(order)_max_BVoigt_st_Vol=$(vol)_dh=$(dh)_gam=$(gamma)_amin=$(amin)_alph=$(alpha_coeff)_ls=$(line_search)$(write_dir_suffix)/",
    iter_mod::Int=10,
    reinit_mod::Int=5,
    solver_method::String="Schur",
    petsc_opts::String="-ksp_converged_reason",
    elast_solver_rtol::Float64=10^-2,
    cond_solver_rtol::Float64=10^-2,
    init_struc::String="")
  mesh_partition = (px,py,pz)
  with_mpi() do distribute
    ranks = distribute(LinearIndices((prod(mesh_partition),)))
    solver_method = Symbol(solver_method)
    if i_am_main(ranks)
      println("--------------------------------------------------")
      println("Parsed args for Max Bₕ s.t., Vol = $vol, dh = $dh:")
      println("       (px,py,pz) => ", mesh_partition)
      println("               dh => ", dh)
      println("              vol => ", vol)
      println("                n => ", n)
      println("            order => ", order)
      println("            gamma => ", gamma)
      println("             amin => ", amin)
      println("      alpha_coeff => ", alpha_coeff)
      println("         lsf_func => ", lsf_func)
      println("        write_dir => ", write_dir)
      println("         iter_mod => ", iter_mod)
      println("       reinit_mod => ", reinit_mod)
      println("    solver_method => ", solver_method)
      println("       petsc_opts => ", petsc_opts)
      println("       init_struc => ", init_struc)
      println("      line_search => ", line_search)
      println("elast_solver_rtol => ", elast_solver_rtol)
      println(" cond_solver_rtol => ", cond_solver_rtol)
      println("--------------------------------------------------")
    end
    lsf_fn(x) = invokelatest(eval(Meta.parse(lsf_func)),x)
    params = (;mesh_partition,dh,vf=vol,n,order,γ=gamma,amin,amax,alpha_coeff,lsf_fn,path=write_dir,
      iter_mod,solver_method,petsc_opts,init_struc,reinit_mod,line_search,xi,xi_reduce_coef,xi_reduce_abs_tol,elast_solver_rtol,cond_solver_rtol)
    GridapPETSc.with(args=split(petsc_opts)) do
      main(ranks,params)
    end
  end
end