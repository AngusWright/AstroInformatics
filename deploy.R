apps <- c(
  "MCMC" = "Metropolis-Hastings MCMC",
  "ShiftFitting" = "Shift Fitting under Heteroskedasticity"
)
output_dir <- "docs"

for (app_dir in names(apps)) {
  shinylive::export(app_dir, output_dir, subdir = app_dir)
}

links <- paste(
  sprintf('      <li><a href="%s/">%s</a></li>', names(apps), unname(apps)),
  collapse = "\n"
)

writeLines(
  c(
    "<!doctype html>",
    '<html lang="en">',
    "  <head>",
    '    <meta charset="utf-8">',
    '    <meta name="viewport" content="width=device-width, initial-scale=1">',
    "    <title>Astroinformatics apps</title>",
    "  </head>",
    "  <body>",
    "    <main>",
    "      <h1>Astroinformatics apps</h1>",
    "      <ul>",
    links,
    "      </ul>",
    "    </main>",
    "  </body>",
    "</html>"
  ),
  file.path(output_dir, "index.html")
)
