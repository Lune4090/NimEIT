import std/[rdstdin, strutils]
import db_connector/db_sqlite
import mesh

proc update_database*(mesh: Mesh, meshName: string) =
  echo "Data is generated, updating database..."
  let
    experimentID = readLineFromStdin("Experiment id: ")
    db = open("data/" & meshName & "/mesh.db", "", "", "")

  echo "Writing database..."

  db.exec(sql"""CREATE TABLE IF NOT EXISTS ElementTable (
                ExperimentID INTEGER,
                ElementID INTEGER,
                σRef FLOAT
            )""")

  db.exec(sql"""CREATE TABLE IF NOT EXISTS VerticeTable (
                ExperimentID INTEGER,
                VerticeID INTEGER,
                J FLOAT,
                V FLOAT
            )""")

  for (i, elem) in mesh.elements.pairs():
    db.exec(sql"INSERT INTO ElementTable (ExperimentID, ElementID, σRef) VALUES (?, ?, ?)",
      $experimentID.parseInt, $i, $elem.σRef)

  for (i, vert) in mesh.vertices.pairs():
    db.exec(sql"INSERT INTO VerticeTable (ExperimentID, VerticeID, J, V) VALUES (?, ?, ?, ?)",
      $experimentID.parseInt, $i, $vert.J, $vert.V)

  echo "Database is updated"

  db.close()
