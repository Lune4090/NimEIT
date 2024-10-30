import std/[math, sugar]
import arraymancer, results
import mesh, calc

proc δσ_over_δV*(jac: Tensor[float], α = 1.0, p = 1.0): Result[Tensor[float], CatchableError] =
  ## https://ieeexplore.ieee.org/document/6971063/
  ## α: Tikhonovの正則化項の係数
  ## p: Newton-Raphson法に基づく正則化行列のスケーリング項
  let
    JtJ = jac.transpose*jac
    Q = (diagonalize(JtJ).value).map(x => x.pow(p))

  return ((JtJ + α^2*Q).pinv * (jac.transpose)).ok()

proc reconstruct_δσ*(mesh: Mesh, coef: Tensor[float]): Result[Tensor[float], CatchableError] =
  var
    δV: seq[float]
  for (i, vert) in mesh.vertices.pairs:
    if i >= mesh.numOuterVertices:
      break
    δV.add(vert.ΔV)

  return (coef*δV.toTensor).ok()

proc compute_jac_2d_tri*(mesh: Mesh, stiffnessMatrix: Tensor[float], stackedLocalStiffnessMatrix: Tensor[float]): Result[Tensor[float], CatchableError] =
  ## pyEITのEITforward.compute_jac()に相当
  if stiffnessMatrix.shape != [len(mesh.vertices), len(mesh.vertices)]:
    return CatchableError(msg: "stiffnessMatrix's shape must be [len(mesh.vertices), len(mesh.vertices)]").err()

  if stackedLocalStiffnessMatrix.shape != [len(mesh.elements), 3, 3]:
    return CatchableError(msg: "stackedLocalStiffnessMatrix's shape must be [len(mesh.elements), 3, 3]").err()

  var
    Vs: seq[float]
    jac = zeros[float]([mesh.numOuterVertices, len(mesh.elements)])
  let
    stiffMatInv = stiffnessMatrix.pinv[0..<mesh.numOuterVertices, _]
  for vertice in mesh.vertices.items:
    Vs.add(vertice.V)

  for (i, element) in mesh.elements.pairs:
    let
      slicedStiffMatInv = concat(stiffMatInv[_, element.idxVertice1], stiffMatInv[_, element.idxVertice2], stiffMatInv[_, element.idxVertice3], axis=1)
      localStiffnessMat = concat(stackedLocalStiffnessMatrix[i, 0, _], stackedLocalStiffnessMatrix[i, 1, _], stackedLocalStiffnessMatrix[i, 2, _], axis=1).reshape(3, 3)
    jac[_, i] = (-slicedStiffMatInv * localStiffnessMat * @[Vs[element.idxVertice1], Vs[element.idxVertice2], Vs[element.idxVertice3]].toTensor)

  return jac.ok()
