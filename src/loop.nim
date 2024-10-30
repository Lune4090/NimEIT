import std/[rdstdin, strutils, sequtils, os, random, tables]
import arraymancer, db_connector/db_sqlite, results
import setting, plotter, backward, mesh, database, toml

proc forward_loop*(preserve_data: bool, meshName: string, settingFileName: string) =

  # Get mesh data and setting data
  let
    meshTomlPath = "data/" & meshName & "/mesh.toml"
    meshParams = mesh_params_from_toml(meshTomlPath).value()
    settingTomlPath = "data/" & meshName & "/" & settingFileName & ".toml"
    (centers, Rs, σRefs) = σRefs_from_toml(settingTomlPath).value()
    (verts, Js) = Js_from_toml(settingTomlPath).value()

  # Generate mesh
  var mesh2d = generate_mesh(meshParams, drawVert = true, drawMesh = true)

  # Elements
  mesh2d.modify_σRef_circle_region(centers, Rs, σRefs)
  mesh2d.modify_J(verts, Js)

  # Get stiffness matrices
  var (stackedLocalStiffnessMat, unitStackedLocalStiffnessMat, stiffness_mat) = get_stiffness_matrices(mesh2d)
  discard stackedLocalStiffnessMat
  discard unitStackedLocalStiffnessMat

  # Forward. Solve KV=J based on Galerkin method and update V, then get the voltage mapping
  var σRef: seq[float]
  for (i, elem) in mesh2d.elements.pairs():
    σRef.add(elem.σRef)

  # 次元はAmpere/Length(ここでスケール反映!)
  var J: seq[float]
  for (i, vert) in mesh2d.vertices.pairs():
    J.add(vert.J)

  let V = solve(stiffness_mat, J.toTensor)
  for (i, vert) in mesh2d.vertices.mpairs():
    vert.V = V[i]

  draw_V(mesh2d, (1000, 1000), ((-meshParams.diameter, -meshParams.diameter), (meshParams.diameter, meshParams.diameter)))

  if preserve_data: update_database(mesh2d, meshName)


proc backward_loop*(meshName: string) =

  # Get mesh data
  let
    meshTomlPath = "data/" & meshName & "/mesh.toml"
    meshParams = mesh_params_from_toml(meshTomlPath).value()

  # Generate mesh
  var mesh2d = generate_mesh(meshParams, drawVert = false, drawMesh = false)

  # Read input
  var inputTomlName = readLineFromStdin("Input .toml file name: ")
  while not fileExists("data/" & meshName & "/" & inputTomlName & ".toml"):
    echo "Input file is not found"
    inputTomlName = readLineFromStdin("Input .toml file name: ")
  var
    experimentIDs0: seq[int]
    experimentIDs1: seq[int]

  (experimentIDs0, experimentIDs1) = experimentIDs_from_toml("data/" & meshName & "/" & inputTomlName & ".toml").value()

  if len(experimentIDs0) != len(experimentIDs1):
    echo "length of experimentID (1st/2nd) is not same, check it again"
    return

  var
    δσs: seq[seq[float]]

  for i in 0..<len(experimentIDs0):
    let
      experimentID0 = experimentIDs0[i]
      experimentID1 = experimentIDs1[i]
    var
      J: seq[float]
      V0: seq[float]
      V1: seq[float]
      σ0: seq[float]
      σ1: seq[float]

    echo "Reading database..."
    let db = open("data/" & meshName & "/mesh.db", "", "", "")

    for row in db.fastRows(sql"SELECT ExperimentID, σRef FROM ElementTable"):
      if row[0].parseInt == experimentID0:
        σ0.add(row[1].parseFloat)
      if row[0].parseInt == experimentID1:
        σ1.add(row[1].parseFloat)

    for row in db.fastRows(sql"SELECT ExperimentID, J, V FROM VerticeTable"):
      if row[0].parseInt == experimentID0:
        J.add(row[1].parseFloat)
        V0.add(row[2].parseFloat)
      if row[0].parseInt == experimentID1:
        V1.add(row[2].parseFloat)

    # ノイズの導入(ここじゃなくてメッシュ本体に直接加算すべきかもしれない、伝導率も同じく)
    echo inputTomlName
    var errors = intentional_error_from_toml("data/" & meshName & "/" & inputTomlName & ".toml").value()
    for error in errors.keys():
      if error == "Vs":
        if errors["Vs"]["type"] == "Gaussian":
          let
            mu = errors["Vs"]["mu"].parseFloat
            sigma = errors["Vs"]["sigma"].parseFloat
          echo "Add gaussian noise to voltages..."
          for V in V0.mitems():
            V = V + gauss(mu = mu, sigma = sigma)
          for V in V1.mitems():
            V = V + gauss(mu = mu, sigma = sigma)

    for error in errors.keys():
      if error == "σs":
        if errors["σs"]["type"] == "Gaussian":
          let
            mu = errors["σs"]["mu"].parseFloat
            sigma = errors["σs"]["sigma"].parseFloat
          echo "Add gaussian noise to voltages..."
          for σ in σ0.mitems():
            σ = σ + gauss(mu = mu, sigma = sigma)
          for σ in σ1.mitems():
            σ = σ + gauss(mu = mu, sigma = sigma)

    db.close()

    for (j, elem) in mesh2d.elements.mpairs():
      elem.σRef = σ0[j]
      elem.Δσ = σ1[j] - σ0[j]
    for (j, vert) in mesh2d.vertices.mpairs():
      vert.J = J[j]
      vert.V = V0[j]
      vert.ΔV = V1[j] - V0[j]

    # Get stiffness matrices
    var (stackedLocalStiffnessMat, unitStackedLocalStiffnessMat, stiffness_mat) = get_stiffness_matrices(mesh2d)
    discard stackedLocalStiffnessMat

    # Backward-1. Calculate jacobian from global / local stiffness matrix and outer node's voltages
    let jac = mesh2d.compute_jac_2d_tri(stiffness_mat, unitStackedLocalStiffnessMat).value

    # Backward-2. Converge RMS based on differential re-construction method with regularization term
    let
      coef = jac.δσ_over_δV().value
      δσ = mesh2d.reconstruct_δσ(coef).value

    var RMS = 0.0
    for (j, elem) in mesh2d.elements.mpairs():
      elem.δσ = δσ[j]
      RMS += sqrt((elem.Δσ - elem.δσ)^2)
    RMS = RMS/len(mesh2d.elements).float

    echo "RMS(" & $i & "): " & $RMS

    # Backward-3. Get the reconstructed image !
    δσs.add(δσ.toSeq1D)

    draw_δσ(mesh2d, (1000, 1000), ((-meshParams.diameter, -meshParams.diameter), (meshParams.diameter, meshParams.diameter)))


  var δσ_mean = repeat(0.0, len(δσs[0]))
  for j in 0..<len(δσs[0]):
    for i in 0..<len(δσs):
      δσ_mean[j] += δσs[i][j]
    δσ_mean[j] = δσ_mean[j]/len(δσs).float

  var RMS = 0.0
  for (i, elem) in mesh2d.elements.mpairs():
    elem.δσ = δσ_mean[i]
    RMS += sqrt((elem.Δσ - elem.δσ)^2)
  RMS = RMS/len(mesh2d.elements).float
  echo "RMS (last): " & $RMS

  draw_Δσ(mesh2d, (1000, 1000), ((-meshParams.diameter, -meshParams.diameter), (meshParams.diameter, meshParams.diameter)))
  draw_δσ(mesh2d, (1000, 1000), ((-meshParams.diameter, -meshParams.diameter), (meshParams.diameter, meshParams.diameter)), title = "δσ_mean(mean of estimated conductivities change)")
