### Weak form for homogenisation problem

# Note 1:
#  - Non-dimensionalisation is based on approach in the paper "Cell-by-cell approximate
#    Schur complement technique in preconditioning of meshfree discretized piezoelectric
#    equations" by Cao and Neytcheva, 2021 (10.1002/nla.2362).
#  - This improves upon the adapted Amlashi approach as it results in better-conditioned RHS.
#  - For all problems here we take the characteristic length scale l=1.
#  - To retrieve tensors in dimensionalised form we use: C=C_nd*c1/A₁₁; e=e_nd*e1/A₁₂; C=C_nd*κ1/A₂₂;
#
# Note 2:
#  - We expect the resulting closure to take arguments (uϕ,φ,dΩ) for the
#    funcitional and (q,uϕ,φ,dΩ) for the derivative.
#  - Inclusion of additional quantities such as Lagrange multipliers require
#    modification of all functionals with hardcoded unpacking of uϕ. These and
#    should be considered seperately in a driver script.
#
# Note 3:
#  - When `sym=false`, we take:
#      |  A    B | * | x | = | b |
#      | -Bᵀ   C |   | y |   | c |
#    This is required for the block Schur preconditioner
#  - When `sym=true`, we take:
#      |  A    B | * | x | = | b |
#      |  Bᵀ  -C |   | y |   |-c |
#    This is required for the TriCG solver

function piezo_hom_weak_form(mat_params::NamedTuple,interp;sym::Bool=true)
  C, e, κ, A₁₁, A₁₂, A₂₂, _, _, _ = mat_params
  I = interp.I
  D = first(size(C))
  s = 1 - 2Int(sym) # -ve if sym is true
  @assert D ∈ (2,3) "Dimension must be 2 or 3"

  if isequal(D,2)
    εᴹ = (SymTensorValue(1.0,0.0,0.0),SymTensorValue(0.0,0.0,1.0),SymTensorValue(0.0,1/2,0))
    Eⁱ = (VectorValue(1.0,0.0,),VectorValue(0.0,1.0))
  else
    εᴹ = (SymTensorValue(1.0,0.0,0.0,0.0,0.0,0.0),
          SymTensorValue(0.0,0.0,0.0,1.0,0.0,0.0),
          SymTensorValue(0.0,0.0,0.0,0.0,0.0,1.0),
          SymTensorValue(0.0,0.0,0.0,0.0,0.5,0.0),
          SymTensorValue(0.0,0.0,0.5,0.0,0.0,0.0),
          SymTensorValue(0.0,0.5,0.0,0.0,0.0,0.0))
    Eⁱ = (VectorValue(1.,0.,0.),VectorValue(0.,1.,0),VectorValue(0.,0.,1.))
  end

  # We take -ve the second row to obtain a symmetric stiffness matrix
  a((u,ϕ),(v,q),φ,dΩ) = ∫((I ∘ φ) * (A₁₁*((C ⊙ ε(u)) ⊙ ε(v))   - A₁₂*((-∇(ϕ) ⋅ e) ⊙ ε(v)) +
                                   s*A₁₂*((e ⋅² ε(u)) ⋅ -∇(q)) + s*A₂₂*((κ ⋅ -∇(ϕ)) ⋅ -∇(q))) )dΩ;

  l_ε = [((v,q),φ,dΩ) -> ∫((I ∘ φ) * (-A₁₁*((C ⊙ εᴹ[i]) ⊙ ε(v)) +
                                    -s*A₁₂*((e ⋅² εᴹ[i]) ⋅ -∇(q))) )dΩ for i = 1:length(εᴹ)]
  l_E = [((v,q),φ,dΩ) -> ∫((I ∘ φ) * ( A₁₂*((Eⁱ[i] ⋅ e) ⊙ ε(v)) +
                                    -s*A₂₂*((κ ⋅ Eⁱ[i]) ⋅ -∇(q))) )dΩ for i = 1:length(Eⁱ)]
  l = [l_ε; l_E]

  return a,l,εᴹ,Eⁱ
end

### Homogenised material tensor functionals

# Note:
#  - In 2D uϕ unpacks as -> uϕ = u1,ϕ1,u2,ϕ2,u3,ϕ3,u4,ϕ4,u5,ϕ5
#    where entries 1,2,3 correspond to εᴹ[1],εᴹ[2],εᴹ[3] and
#    4,5 correspond to Eⁱ[1],Eⁱ[2].
#  - In 3D uϕ unpacks as -> uϕ = u1,ϕ1,u2,ϕ2,u3,ϕ3,u4,ϕ4,u5,ϕ5,u6,ϕ6,u7,ϕ7,u8,ϕ8,u9,ϕ9
#    where entries 1,2,3,4,5,6 correspond to εᴹ[1],εᴹ[2],εᴹ[3],εᴹ[4],εᴹ[5],εᴹ[6] and
#    7,8,9 correspond to Eⁱ[1],Eⁱ[2],Eⁱ[3].

function homogenised_mat_tensor_functionals(mat_params::NamedTuple,interp,εᴹ,Eⁱ)
  C, e, κ, A₁₁, A₁₂, A₂₂, _, _, _ = mat_params
  I, DH = interp.I, interp.DH
  N = 3*(first(size(C))-1) # <- Number of unique strain solutions

  function Cᴴ(r,s,uϕ,φ,dΩ)
    u_s = uϕ[2s-1]; ϕ_s = uϕ[2s]
    ∫((I ∘ φ) * (A₁₁*((C ⊙ (ε(u_s) + εᴹ[s])) ⊙ εᴹ[r]) - A₁₂*((-∇(ϕ_s) ⋅ e) ⊙ εᴹ[r])))dΩ;
  end
  function DCᴴ(r,s,q,uϕ,φ,dΩ)
    u_r = uϕ[2r-1]; ϕ_r = uϕ[2r]
    u_s = uϕ[2s-1]; ϕ_s = uϕ[2s]
    ∫(- q * (
      A₁₁*((C ⊙ (ε(u_s) + εᴹ[s])) ⊙ (ε(u_r) + εᴹ[r])) -
      A₁₂*((-∇(ϕ_s) ⋅ e) ⊙ (ε(u_r) + εᴹ[r])) -
      A₁₂*((e ⋅² (ε(u_s) + εᴹ[s])) ⋅ (-∇(ϕ_r))) -
      A₂₂*((κ ⋅ (-∇(ϕ_s))) ⋅ (-∇(ϕ_r)))
      ) * (DH ∘ φ) * (norm ∘ ∇(φ))
    )dΩ;
  end

  function eᴴ(i,s,uϕ,φ,dΩ)
    u_s = uϕ[2s-1]; ϕ_s = uϕ[2s]
    ∫( (I ∘ φ) * (A₁₂*((e ⋅² (ε(u_s) + εᴹ[s])) ⋅ Eⁱ[i]) + A₂₂*((κ ⋅ (-∇(ϕ_s))) ⋅ Eⁱ[i])))dΩ;
  end
  function Deᴴ(i,s,q,uϕ,φ,dΩ)
    u_i = uϕ[2N+2i-1]; ϕ_i = uϕ[2N+2i]
    u_s = uϕ[2s-1]; ϕ_s = uϕ[2s]
    ∫(-  q * (
      - A₁₁*((C ⊙ (ε(u_s) + εᴹ[s])) ⊙ (ε(u_i))) +
        A₁₂*((-∇(ϕ_s) ⋅ e) ⊙ (ε(u_i))) +
        A₁₂*((e ⋅² (ε(u_s) + εᴹ[s])) ⋅ (-∇(ϕ_i) + Eⁱ[i])) +
        A₂₂*((κ ⋅ (-∇(ϕ_s))) ⋅ (-∇(ϕ_i) + Eⁱ[i]))
      ) * (DH ∘ φ) * (norm ∘ ∇(φ))
    )dΩ;
  end

  function κᴴ(i,j,uϕ,φ,dΩ)
    u_j = uϕ[2N+2j-1]; ϕ_j = uϕ[2N+2j]
    ∫((I ∘ φ) * (A₁₂*((e ⋅² ε(u_j)) ⋅ Eⁱ[i]) + A₂₂*((κ ⋅ (-∇(ϕ_j) + Eⁱ[j])) ⋅ Eⁱ[i])))dΩ;
  end
  function Dκᴴ(i,j,q,uϕ,φ,dΩ)
    u_i = uϕ[2N+2i-1]; ϕ_i = uϕ[2N+2i]
    u_j = uϕ[2N+2j-1]; ϕ_j = uϕ[2N+2j]
    ∫(- q * (
      - A₁₁*((C ⊙ (ε(u_j))) ⊙ (ε(u_i))) +
        A₁₂*(((-∇(ϕ_j) + Eⁱ[j]) ⋅ e) ⊙ (ε(u_i))) +
        A₁₂*((e ⋅² (ε(u_j))) ⋅ (-∇(ϕ_i) + Eⁱ[i])) +
        A₂₂*((κ ⋅ (-∇(ϕ_j) + Eⁱ[j])) ⋅ (-∇(ϕ_i) + Eⁱ[i]))
      ) * (DH ∘ φ) * (norm ∘ ∇(φ))
    )dΩ;
  end

  return Cᴴ,eᴴ,κᴴ,DCᴴ,Deᴴ,Dκᴴ
end

#### Optimisation functionals
using ForwardDiff
get_value(x) = x
get_value(x::ForwardDiff.Dual) = x.value

function bulk_modulus_voigt_2d(Cᴴ,eᴴ,κᴴ,DCᴴ,Deᴴ,Dκᴴ)
  Bᴴ(uϕ,φ,dΩ) = 1/4*(Cᴴ(1,1,uϕ,φ,dΩ)+Cᴴ(2,2,uϕ,φ,dΩ)+2*Cᴴ(1,2,uϕ,φ,dΩ))
  DBᴴ(q,uϕ,φ,dΩ) = 1/4*(DCᴴ(1,1,q,uϕ,φ,dΩ)+DCᴴ(2,2,q,uϕ,φ,dΩ)+2*DCᴴ(1,2,q,uϕ,φ,dΩ))

  return Bᴴ,DBᴴ
end

function bulk_modulus_voigt_3d(Cᴴ,eᴴ,κᴴ,DCᴴ,Deᴴ,Dκᴴ)
  Bᴴ(uϕ,φ,dΩ) = 1/9*(Cᴴ(1,1,uϕ,φ,dΩ)+Cᴴ(2,2,uϕ,φ,dΩ)+Cᴴ(3,3,uϕ,φ,dΩ)+
    2*(Cᴴ(1,2,uϕ,φ,dΩ)+Cᴴ(1,3,uϕ,φ,dΩ)+Cᴴ(2,3,uϕ,φ,dΩ)))
  DBᴴ(q,uϕ,φ,dΩ) = 1/9*(DCᴴ(1,1,q,uϕ,φ,dΩ)+DCᴴ(2,2,q,uϕ,φ,dΩ)+DCᴴ(3,3,q,uϕ,φ,dΩ)+
    2*(DCᴴ(1,2,q,uϕ,φ,dΩ)+DCᴴ(1,3,q,uϕ,φ,dΩ)+DCᴴ(2,3,q,uϕ,φ,dΩ)))

  return Bᴴ,DBᴴ
end

function bulk_modulus_2d(Cᴴ,eᴴ,κᴴ,DCᴴ,Deᴴ,Dκᴴ)
  idxs = [CartesianIndex(j,i) for i = 1:3 for j = i:3]
  _Cᴴ = zeros(Float64,3,3)

  ## Bulk modulus
  function Bᴴ(uϕ,φ,dΩ)
    _Cᴴ11 = sum(Cᴴ(idxs[1][1],idxs[1][2],uϕ,φ,dΩ))
    _Cᴴ[idxs[1]] = get_value(_Cᴴ11)
    for idx in idxs[2:end]
      _Cᴴ[idx] = get_value(sum(Cᴴ(idx[1],idx[2],uϕ,φ,dΩ)));
    end
    Sᴴ = inv(Symmetric(_Cᴴ,:L))

    return (1/sum(Sᴴ[i,j] for i = 1:2, j = 1:2))*∫(1)dΩ # assume volume is 1!!!
  end

  function DBᴴ(q,uϕ,φ,dΩ)
    Sᴴ = inv(Symmetric(_Cᴴ,:L))
    Cᴴ′ = map(idx->DCᴴ(idx[1],idx[2],q,uϕ,φ,dΩ),CartesianIndices((1:3,1:3)))

    Bᴴ = 1/sum(Sᴴ[i,j] for i = 1:2, j = 1:2)

    # Bᴴ′ = 1/Sᴴᵢᵢⱼⱼ²(SᴴᵢᵢₖₖCᴴ′ₖₖₗₗSᴴₗₗₘₘ)
    return ((Sᴴ[1,1]+Sᴴ[2,1])^2*Bᴴ^2)*Cᴴ′[1,1]+(2*(Sᴴ[1,1]+Sᴴ[2,1])*(Sᴴ[2,1]+Sᴴ[2,2])*Bᴴ^2)*Cᴴ′[2,1]+
      ((Sᴴ[2,1]+Sᴴ[2,2])^2*Bᴴ^2)*Cᴴ′[2,2] + (2*(Sᴴ[1,1]+Sᴴ[2,1])*(Sᴴ[3,1]+Sᴴ[3,2])*Bᴴ^2)*Cᴴ′[3,1]+
      (2*(Sᴴ[2,1]+Sᴴ[2,2])*(Sᴴ[3,1]+Sᴴ[3,2])*Bᴴ^2)*Cᴴ′[3,2]+((Sᴴ[3,1]+Sᴴ[3,2])^2*Bᴴ^2)*Cᴴ′[3,3]
  end

  return Bᴴ,DBᴴ
end

function bulk_modulus_3d(Cᴴ,eᴴ,κᴴ,DCᴴ,Deᴴ,Dκᴴ)
  idxs = [CartesianIndex(j,i) for i = 1:6 for j = i:6]
  _Cᴴ = zeros(Float64,6,6)

  function Bᴴ(uϕ,φ,dΩ)
    _Cᴴ11 = sum(Cᴴ(idxs[1][1],idxs[1][2],uϕ,φ,dΩ))
    _Cᴴ[idxs[1]] = get_value(_Cᴴ11)
    for idx in idxs[2:end]
      _Cᴴ[idx] = get_value(sum(Cᴴ(idx[1],idx[2],uϕ,φ,dΩ)));
    end
    Sᴴ = inv(Symmetric(_Cᴴ,:L))

    return (1/sum(Sᴴ[i,j] for i = 1:3, j = 1:3))*∫(1)dΩ
  end

  function DBᴴ(q,uϕ,φ,dΩ)
    Sᴴ = inv(Symmetric(_Cᴴ,:L))
    Cᴴ′ = map(idx->DCᴴ(idx[1],idx[2],q,uϕ,φ,dΩ),CartesianIndices((1:6,1:6)))

    Bᴴ = 1/sum(Sᴴ[i,j] for i = 1:3, j = 1:3)

    # Bᴴ′ = 1/Sᴴᵢᵢⱼⱼ²(SᴴᵢᵢₖₖCᴴ′ₖₖₗₗSᴴₗₗₘₘ)
    return (Bᴴ^2*(Sᴴ[1,1]+Sᴴ[2,1]+Sᴴ[3,1])^2)*Cᴴ′[1,1]+(2*Bᴴ^2*(Sᴴ[1,1]+Sᴴ[2,1]+Sᴴ[3,1])*(Sᴴ[2,1]+Sᴴ[2,2]+Sᴴ[3,2]))*Cᴴ′[2,1]+(Bᴴ^2*(Sᴴ[2,1]+Sᴴ[2,2]+Sᴴ[3,2])^2)*Cᴴ′[2,2]+(2*Bᴴ^2*(Sᴴ[1,1]+Sᴴ[2,1]+Sᴴ[3,1])*(Sᴴ[3,1]+Sᴴ[3,2]+Sᴴ[3,3]))*Cᴴ′[3,1]+(2*Bᴴ^2*(Sᴴ[2,1]+Sᴴ[2,2]+Sᴴ[3,2])*(Sᴴ[3,1]+Sᴴ[3,2]+Sᴴ[3,3]))*Cᴴ′[3,2]+(Bᴴ^2*(Sᴴ[3,1]+Sᴴ[3,2]+Sᴴ[3,3])^2)*Cᴴ′[3,3]+(2*Bᴴ^2*(Sᴴ[1,1]+Sᴴ[2,1]+Sᴴ[3,1])*(Sᴴ[4,1]+Sᴴ[4,2]+Sᴴ[4,3]))*Cᴴ′[4,1]+(2*Bᴴ^2*(Sᴴ[2,1]+Sᴴ[2,2]+Sᴴ[3,2])*(Sᴴ[4,1]+Sᴴ[4,2]+Sᴴ[4,3]))*Cᴴ′[4,2]+(2*Bᴴ^2*(Sᴴ[3,1]+Sᴴ[3,2]+Sᴴ[3,3])*(Sᴴ[4,1]+Sᴴ[4,2]+Sᴴ[4,3]))*Cᴴ′[4,3]+(Bᴴ^2*(Sᴴ[4,1]+Sᴴ[4,2]+Sᴴ[4,3])^2)*Cᴴ′[4,4]+(2*Bᴴ^2*(Sᴴ[1,1]+Sᴴ[2,1]+Sᴴ[3,1])*(Sᴴ[5,1]+Sᴴ[5,2]+Sᴴ[5,3]))*Cᴴ′[5,1]+(2*Bᴴ^2*(Sᴴ[2,1]+Sᴴ[2,2]+Sᴴ[3,2])*(Sᴴ[5,1]+Sᴴ[5,2]+Sᴴ[5,3]))*Cᴴ′[5,2]+(2*Bᴴ^2*(Sᴴ[3,1]+Sᴴ[3,2]+Sᴴ[3,3])*(Sᴴ[5,1]+Sᴴ[5,2]+Sᴴ[5,3]))*Cᴴ′[5,3]+(2*Bᴴ^2*(Sᴴ[4,1]+Sᴴ[4,2]+Sᴴ[4,3])*(Sᴴ[5,1]+Sᴴ[5,2]+Sᴴ[5,3]))*Cᴴ′[5,4]+(Bᴴ^2*(Sᴴ[5,1]+Sᴴ[5,2]+Sᴴ[5,3])^2)*Cᴴ′[5,5]+(2*Bᴴ^2*(Sᴴ[1,1]+Sᴴ[2,1]+Sᴴ[3,1])*(Sᴴ[6,1]+Sᴴ[6,2]+Sᴴ[6,3]))*Cᴴ′[6,1]+(2*Bᴴ^2*(Sᴴ[2,1]+Sᴴ[2,2]+Sᴴ[3,2])*(Sᴴ[6,1]+Sᴴ[6,2]+Sᴴ[6,3]))*Cᴴ′[6,2]+(2*Bᴴ^2*(Sᴴ[3,1]+Sᴴ[3,2]+Sᴴ[3,3])*(Sᴴ[6,1]+Sᴴ[6,2]+Sᴴ[6,3]))*Cᴴ′[6,3]+(2*Bᴴ^2*(Sᴴ[4,1]+Sᴴ[4,2]+Sᴴ[4,3])*(Sᴴ[6,1]+Sᴴ[6,2]+Sᴴ[6,3]))*Cᴴ′[6,4]+(2*Bᴴ^2*(Sᴴ[5,1]+Sᴴ[5,2]+Sᴴ[5,3])*(Sᴴ[6,1]+Sᴴ[6,2]+Sᴴ[6,3]))*Cᴴ′[6,5]+(Bᴴ^2*(Sᴴ[6,1]+Sᴴ[6,2]+Sᴴ[6,3])^2)*Cᴴ′[6,6]
  end

  return Bᴴ,DBᴴ
end

## Hydrostatic coupling

# This can be replaced with an @generated function to compute all of these quantities.
function hydrostatic_coupling_2d(Cᴴ,eᴴ,κᴴ,DCᴴ,Deᴴ,Dκᴴ)
  _Cᴴ = zeros(Float64,3,3)
  _eᴴ = zeros(Float64,2,3)
  idxs = [CartesianIndex(j,i) for i = 1:3 for j = i:3]

  function dᴴ(uϕ,φ,dΩ)
    _Cᴴ11 = sum(Cᴴ(idxs[1][1],idxs[1][2],uϕ,φ,dΩ))
    _Cᴴ[idxs[1]] = get_value(_Cᴴ11)
    for idx in idxs[2:end]
      _Cᴴ[idx] = get_value(sum(Cᴴ(idx[1],idx[2],uϕ,φ,dΩ)));
    end
    Sᴴ = inv(Symmetric(_Cᴴ,:L))
    eᴴ21 = eᴴ(2,1,uϕ,φ,dΩ);
    eᴴ22 = eᴴ(2,2,uϕ,φ,dΩ);
    eᴴ23 = eᴴ(2,3,uϕ,φ,dΩ);

    # dᴴ = eᴴ*inv(Cᴴ); dₕ = dᴴ[2,1] + dᴴ[2,2]
    return (Sᴴ[1,1]+Sᴴ[2,1])*eᴴ21+(Sᴴ[2,1]+Sᴴ[2,2])*eᴴ22+(Sᴴ[3,1]+Sᴴ[3,2])*eᴴ23
  end
  function Ddᴴ(q,uϕ,φ,dΩ)
    ## Assume we're always recalculating dᴴ when finding Ddᴴ # IGNORE
    _Cᴴ11 = sum(Cᴴ(idxs[1][1],idxs[1][2],uϕ,φ,dΩ))
    _Cᴴ[idxs[1]] = get_value(_Cᴴ11)
    for idx in idxs[2:end]
      _Cᴴ[idx] = get_value(sum(Cᴴ(idx[1],idx[2],uϕ,φ,dΩ)));
    end
    ## Other entries don't matter in computation
    _eᴴ[2,1] = get_value(sum(eᴴ(2,1,uϕ,φ,dΩ)));
    _eᴴ[2,2] = get_value(sum(eᴴ(2,2,uϕ,φ,dΩ)));
    _eᴴ[2,3] = get_value(sum(eᴴ(2,3,uϕ,φ,dΩ)));

    Sᴴ = inv(Symmetric(_Cᴴ,:L))
    eS = _eᴴ*Sᴴ

    Cᴴ′ = map(idx->DCᴴ(idx[1],idx[2],q,uϕ,φ,dΩ),CartesianIndices((1:3,1:3)))
    eᴴ′ = map(idx->Deᴴ(idx[1],idx[2],q,uϕ,φ,dΩ),CartesianIndices((1:2,1:3)))

    ## Short but not optimised # Requires S/MMatrix above for _Cᴴ and _eᴴ
    # Cᴴ′ = SMatrix{3,3}(map(idx->DCᴴ(idx[1],idx[2],q,uϕ,φ,dΩ),CartesianIndices((1:3,1:3))))
    # eᴴ′ = SMatrix{2,3}(map(idx->Deᴴ(idx[1],idx[2],q,uϕ,φ,dΩ),CartesianIndices((1:2,1:3))))
    # Ddᴴ = eᴴ′*sym_Sᴴ - eS*Cᴴ′*sym_Sᴴ
    # return Ddᴴ[2,1] + Ddᴴ[2,2]
    ## Long but shit
    return (-eS[2,1]*(Sᴴ[1,1]+Sᴴ[2,1]))*Cᴴ′[1,1]+(-eS[2,2]*(Sᴴ[1,1]+Sᴴ[2,1])-eS[2,1]*(Sᴴ[2,1]+Sᴴ[2,2]))*Cᴴ′[2,1]-
      (eS[2,2]*(Sᴴ[2,1]+Sᴴ[2,2]))*Cᴴ′[2,2]+(-eS[2,3]*(Sᴴ[1,1]+Sᴴ[2,1])-eS[2,1]*(Sᴴ[3,1]+Sᴴ[3,2]))*Cᴴ′[3,1]+
      (-eS[2,3]*(Sᴴ[2,1]+Sᴴ[2,2])-eS[2,2]*(Sᴴ[3,1]+Sᴴ[3,2]))*Cᴴ′[3,2]-eS[2,3]*(Sᴴ[3,1]+Sᴴ[3,2])*Cᴴ′[3,3]+
      (Sᴴ[1,1]+Sᴴ[2,1])*eᴴ′[2,1]+(Sᴴ[2,1]+Sᴴ[2,2])*eᴴ′[2,2]+(Sᴴ[3,1]+Sᴴ[3,2])*eᴴ′[2,3]
  end

  return dᴴ,Ddᴴ
end

function hydrostatic_coupling_3d(Cᴴ,eᴴ,κᴴ,DCᴴ,Deᴴ,Dκᴴ)
  idxs = [CartesianIndex(j,i) for i = 1:6 for j = i:6]
  _Cᴴ = zeros(Float64,6,6)
  _eᴴ = zeros(Float64,3,6)

  function dᴴ(uϕ,φ,dΩ)
    _Cᴴ11 = sum(Cᴴ(idxs[1][1],idxs[1][2],uϕ,φ,dΩ))
    _Cᴴ[idxs[1]] = get_value(_Cᴴ11)
    for idx in idxs[2:end]
      _Cᴴ[idx] = get_value(sum(Cᴴ(idx[1],idx[2],uϕ,φ,dΩ)));
    end
    Sᴴ = inv(Symmetric(_Cᴴ,:L))
    eᴴ31 = eᴴ(3,1,uϕ,φ,dΩ);
    eᴴ32 = eᴴ(3,2,uϕ,φ,dΩ);
    eᴴ33 = eᴴ(3,3,uϕ,φ,dΩ);
    eᴴ34 = eᴴ(3,4,uϕ,φ,dΩ);
    eᴴ35 = eᴴ(3,5,uϕ,φ,dΩ);
    eᴴ36 = eᴴ(3,6,uϕ,φ,dΩ);

    # dᴴ = eᴴ*inv(Cᴴ); dₕ = dᴴ[3,1] + dᴴ[3,2] + dᴴ[3,3]
    return eᴴ31*(Sᴴ[1,1]+Sᴴ[2,1]+Sᴴ[3,1])+eᴴ32*(Sᴴ[2,1]+Sᴴ[2,2]+Sᴴ[3,2])+
      eᴴ33*(Sᴴ[3,1]+Sᴴ[3,2]+Sᴴ[3,3])+eᴴ34*(Sᴴ[4,1]+Sᴴ[4,2]+Sᴴ[4,3])+
      eᴴ35*(Sᴴ[5,1]+Sᴴ[5,2]+Sᴴ[5,3])+eᴴ36*(Sᴴ[6,1]+Sᴴ[6,2]+Sᴴ[6,3])
  end

  function Ddᴴ(q,uϕ,φ,dΩ)
    _Cᴴ11 = sum(Cᴴ(idxs[1][1],idxs[1][2],uϕ,φ,dΩ))
    _Cᴴ[idxs[1]] = get_value(_Cᴴ11)
    for idx in idxs[2:end]
      _Cᴴ[idx] = get_value(sum(Cᴴ(idx[1],idx[2],uϕ,φ,dΩ)));
    end
    _eᴴ[3,1] = get_value(sum(eᴴ(3,1,uϕ,φ,dΩ)));
    _eᴴ[3,2] = get_value(sum(eᴴ(3,2,uϕ,φ,dΩ)));
    _eᴴ[3,3] = get_value(sum(eᴴ(3,3,uϕ,φ,dΩ)));
    _eᴴ[3,4] = get_value(sum(eᴴ(3,4,uϕ,φ,dΩ)));
    _eᴴ[3,5] = get_value(sum(eᴴ(3,5,uϕ,φ,dΩ)));
    _eᴴ[3,6] = get_value(sum(eᴴ(3,6,uϕ,φ,dΩ)));

    Sᴴ = inv(Symmetric(_Cᴴ,:L))
    eS = _eᴴ*Sᴴ

    Cᴴ′ = map(idx->DCᴴ(idx[1],idx[2],q,uϕ,φ,dΩ),CartesianIndices((1:6,1:6)))
    eᴴ′ = map(idx->Deᴴ(idx[1],idx[2],q,uϕ,φ,dΩ),CartesianIndices((1:3,1:6)))

    # dᴴ' = eᴴ′Sᴴ - eᴴSᴴCᴴ′Sᴴ; dₕ' = dᴴ'[3,1] + dᴴ'[3,2] + dᴴ'[3,3]
    return (-eS[3,1]*(Sᴴ[1,1]+Sᴴ[2,1]+Sᴴ[3,1]))*Cᴴ′[1,1]+(-eS[3,2]*(Sᴴ[1,1]+Sᴴ[2,1]+Sᴴ[3,1])-eS[3,1]*(Sᴴ[2,1]+Sᴴ[2,2]+Sᴴ[3,2]))*Cᴴ′[2,1]-
      (eS[3,2]*(Sᴴ[2,1]+Sᴴ[2,2]+Sᴴ[3,2]))*Cᴴ′[2,2]+(-eS[3,3]*(Sᴴ[1,1]+Sᴴ[2,1]+Sᴴ[3,1])-eS[3,1]*(Sᴴ[3,1]+Sᴴ[3,2]+Sᴴ[3,3]))*Cᴴ′[3,1]+
      (-eS[3,3]*(Sᴴ[2,1]+Sᴴ[2,2]+Sᴴ[3,2])-eS[3,2]*(Sᴴ[3,1]+Sᴴ[3,2]+Sᴴ[3,3]))*Cᴴ′[3,2]-(eS[3,3]*(Sᴴ[3,1]+Sᴴ[3,2]+Sᴴ[3,3]))*Cᴴ′[3,3]+
      (-eS[3,4]*(Sᴴ[1,1]+Sᴴ[2,1]+Sᴴ[3,1])-eS[3,1]*(Sᴴ[4,1]+Sᴴ[4,2]+Sᴴ[4,3]))*Cᴴ′[4,1]+
      (-eS[3,4]*(Sᴴ[2,1]+Sᴴ[2,2]+Sᴴ[3,2])-eS[3,2]*(Sᴴ[4,1]+Sᴴ[4,2]+Sᴴ[4,3]))*Cᴴ′[4,2]+
      (-eS[3,4]*(Sᴴ[3,1]+Sᴴ[3,2]+Sᴴ[3,3])-eS[3,3]*(Sᴴ[4,1]+Sᴴ[4,2]+Sᴴ[4,3]))*Cᴴ′[4,3]-(eS[3,4]*(Sᴴ[4,1]+Sᴴ[4,2]+Sᴴ[4,3]))*Cᴴ′[4,4]+
      (-eS[3,5]*(Sᴴ[1,1]+Sᴴ[2,1]+Sᴴ[3,1])-eS[3,1]*(Sᴴ[5,1]+Sᴴ[5,2]+Sᴴ[5,3]))*Cᴴ′[5,1]+
      (-eS[3,5]*(Sᴴ[2,1]+Sᴴ[2,2]+Sᴴ[3,2])-eS[3,2]*(Sᴴ[5,1]+Sᴴ[5,2]+Sᴴ[5,3]))*Cᴴ′[5,2]+
      (-eS[3,5]*(Sᴴ[3,1]+Sᴴ[3,2]+Sᴴ[3,3])-eS[3,3]*(Sᴴ[5,1]+Sᴴ[5,2]+Sᴴ[5,3]))*Cᴴ′[5,3]+
      (-eS[3,5]*(Sᴴ[4,1]+Sᴴ[4,2]+Sᴴ[4,3])-eS[3,4]*(Sᴴ[5,1]+Sᴴ[5,2]+Sᴴ[5,3]))*Cᴴ′[5,4]-(eS[3,5]*(Sᴴ[5,1]+Sᴴ[5,2]+Sᴴ[5,3]))*Cᴴ′[5,5]+
      (-eS[3,6]*(Sᴴ[1,1]+Sᴴ[2,1]+Sᴴ[3,1])-eS[3,1]*(Sᴴ[6,1]+Sᴴ[6,2]+Sᴴ[6,3]))*Cᴴ′[6,1]+
      (-eS[3,6]*(Sᴴ[2,1]+Sᴴ[2,2]+Sᴴ[3,2])-eS[3,2]*(Sᴴ[6,1]+Sᴴ[6,2]+Sᴴ[6,3]))*Cᴴ′[6,2]+
      (-eS[3,6]*(Sᴴ[3,1]+Sᴴ[3,2]+Sᴴ[3,3])-eS[3,3]*(Sᴴ[6,1]+Sᴴ[6,2]+Sᴴ[6,3]))*Cᴴ′[6,3]+
      (-eS[3,6]*(Sᴴ[4,1]+Sᴴ[4,2]+Sᴴ[4,3])-eS[3,4]*(Sᴴ[6,1]+Sᴴ[6,2]+Sᴴ[6,3]))*Cᴴ′[6,4]+
      (-eS[3,6]*(Sᴴ[5,1]+Sᴴ[5,2]+Sᴴ[5,3])-eS[3,5]*(Sᴴ[6,1]+Sᴴ[6,2]+Sᴴ[6,3]))*Cᴴ′[6,5]-(eS[3,6]*(Sᴴ[6,1]+Sᴴ[6,2]+Sᴴ[6,3]))*Cᴴ′[6,6]+
      (Sᴴ[1,1]+Sᴴ[2,1]+Sᴴ[3,1])*eᴴ′[3,1]+(Sᴴ[2,1]+Sᴴ[2,2]+Sᴴ[3,2])*eᴴ′[3,2]+(Sᴴ[3,1]+Sᴴ[3,2]+Sᴴ[3,3])*eᴴ′[3,3]+
      (Sᴴ[4,1]+Sᴴ[4,2]+Sᴴ[4,3])*eᴴ′[3,4]+(Sᴴ[5,1]+Sᴴ[5,2]+Sᴴ[5,3])*eᴴ′[3,5]+(Sᴴ[6,1]+Sᴴ[6,2]+Sᴴ[6,3])*eᴴ′[3,6]
  end

  return dᴴ,Ddᴴ
end