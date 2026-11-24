node <- Sys.which("node")
if (!nzchar(node)) stop("Install Node.js 22+ and run `node --test TeachingTools/test_export.js`.", call. = FALSE)
status <- system2(node, c("--test", "TeachingTools/test_export.js"))
if (status != 0L) stop("HTML export tests failed with exit status ", status, call. = FALSE)
