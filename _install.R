dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE, showWarnings = FALSE)
install.packages("deSolve", repos = "https://cran.r-project.org", lib = Sys.getenv("R_LIBS_USER"), quiet = TRUE)
cat("deSolve installed successfully\n")
