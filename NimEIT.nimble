# Package

version       = "0.1.0"
author        = "lune4090"
description   = "A new awesome nimble package"
license       = "MIT"
srcDir        = "src"
bin           = @["NimEIT"]

# Dependencies

requires "nim >= 2.0.4", "arraymancer", "parsetoml", "db_connector", "plotly", "results", "serial"
