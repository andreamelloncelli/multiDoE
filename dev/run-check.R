devtools::check(
    pkg = ".",
    document = TRUE,
    # build_args = NULL,
    # manual = FALSE,
    # cran = TRUE,
    # remote = FALSE,
    # force_suggests = FALSE,
    run_dont_test = TRUE,
    args = c(
        "--timings",
        "--no-examples"
        # "--no-tests"
    ),
    # env_vars = c(NOT_CRAN = "true"),
    # quiet = FALSE,
    # check_dir = NULL,
    # cleanup = deprecated(),
    vignettes = FALSE,
    error_on = "note",  # c("never", "error", "warning", "note")

)
