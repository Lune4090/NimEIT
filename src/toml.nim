import tables
import results, parsetoml
import setting

proc mesh_params_from_toml*(path: string): Result[MeshParams, CatchableError] =
  let table = parseFile(path)
  if not table["params"].hasKey("numElectrodes"):
    return CatchableError(msg: ".toml format is invalid, numElectrodes is not found.").err()
  if not table["params"].hasKey("diameter"):
    return CatchableError(msg: ".toml format is invalid, diameter is not found.").err()
  if not table["params"].hasKey("numElectrodes"):
    return CatchableError(msg: ".toml format is invalid, numsInnerVertices is not found.").err()
  if not table["params"].hasKey("numElectrodes"):
    return CatchableError(msg: ".toml format is invalid, diameters is not found.").err()
  var
    numElectrodes = table["params"]["numElectrodes"].getInt
    diameter = table["params"]["diameter"].getFloat
    numsInnerVertices: seq[int]
    diameters: seq[float]

  for i in 0..<len(table["params"]["numsInnerVertices"]):
    numsInnerVertices.add(table["params"]["numsInnerVertices"][i].getInt)
  for i in 0..<len(table["params"]["diameters"]):
    diameters.add(table["params"]["diameters"][i].getFloat)

  return MeshParams(numElectrodes: numElectrodes, diameter: diameter, numsInnerVertices: numsInnerVertices, diameters: diameters).ok()

proc σRefs_from_toml*(path: string): Result[(seq[(float, float)], seq[float], seq[float]), CatchableError] =
  let table = parseFile(path)
  if not table["sigmas"].hasKey("centers"):
    return CatchableError(msg: ".toml format is invalid, centers is not found.").err()
  if not table["sigmas"].hasKey("Rs"):
    return CatchableError(msg: ".toml format is invalid, Rs is not found.").err()
  if not table["sigmas"].hasKey("sigmaRefs"):
    return CatchableError(msg: ".toml format is invalid, sigmaRefs is not found.").err()
  var
    centers: seq[(float, float)]
    Rs: seq[float]
    σRefs: seq[float]

  for i in 0..<len(table["sigmas"]["centers"]):
    centers.add((table["sigmas"]["centers"][i][0].getFloat, table["sigmas"]["centers"][i][1].getFloat))
  for i in 0..<len(table["sigmas"]["Rs"]):
    Rs.add(table["sigmas"]["Rs"][i].getFloat)
  for i in 0..<len(table["sigmas"]["sigmaRefs"]):
    σRefs.add(table["sigmas"]["sigmaRefs"][i].getFloat)

  return (centers, Rs, σRefs).ok()

proc Js_from_toml*(path: string): Result[(seq[int], seq[float]), CatchableError] =
  let table = parseFile(path)
  if not table["Js"].hasKey("verts"):
    return CatchableError(msg: ".toml format is invalid, verts is not found.").err()
  if not table["Js"].hasKey("Js"):
    return CatchableError(msg: ".toml format is invalid, Js is not found.").err()

  var
    verts: seq[int]
    Js: seq[float]

  for i in 0..<len(table["Js"]["verts"]):
    verts.add(table["Js"]["verts"][i].getInt)
  for i in 0..<len(table["Js"]["Js"]):
    Js.add(table["Js"]["Js"][i].getFloat)

  return (verts, Js).ok()

proc experimentIDs_from_toml*(path: string): Result[(seq[int], seq[int]), CatchableError] =
  let table = parseFile(path)
  if not table["input"].hasKey("1stExperimentIDs"):
    return CatchableError(msg: ".toml format is invalid, 1stExperimentIDs is not found.").err()
  if not table["input"].hasKey("2ndExperimentIDs"):
    return CatchableError(msg: ".toml format is invalid, 2ndExperimentIDs is not found.").err()

  var
    experimentIDs0: seq[int]
    experimentIDs1: seq[int]

  for i in 0..<len(table["input"]["1stExperimentIDs"]):
    experimentIDs0.add(table["input"]["1stExperimentIDs"][i].getInt)
  for i in 0..<len(table["input"]["2ndExperimentIDs"]):
    experimentIDs1.add(table["input"]["2ndExperimentIDs"][i].getInt)

  return (experimentIDs0, experimentIDs1).ok()

proc intentional_error_from_toml*(path: string): Result[(Table[string, Table[string, string]]), CatchableError] =
  let table = parseFile(path)
  if table.hasKey("error"):
    var
      errors: Table[string, Table[string, string]]

    # エラーを導入する
    if table["error"].hasKey("Vs"):
      echo "Inducing a posterior errors..."
      # とりあえずガウシアンノイズだけ実装
      if table["error"]["Vs"].hasKey("Gaussian"):
        errors["Vs"] = {
          "type": "Gaussian",
          "mu": $table["error"]["Vs"]["Gaussian"]["mu"].getFloat,
          "sigma": $table["error"]["Vs"]["Gaussian"]["sigma"].getFloat,
        }.toTable

    return errors.ok()
  else:
    var
      errors: Table[string, Table[string, string]]
    errors["Vs"] = {
      "type": "Gaussian",
      "mu": $0.0,
      "sigma": $0.0,
    }.toTable
    return errors.ok()
