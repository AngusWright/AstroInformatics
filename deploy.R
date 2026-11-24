node <- Sys.which("node")
if (!nzchar(node)) stop("Install Node.js 22+ and run `node deploy.js`.", call. = FALSE)
status <- system2(node, "deploy.js")
if (status != 0L) stop("HTML site build failed with exit status ", status, call. = FALSE)
