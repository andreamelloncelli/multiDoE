devtools::check(
    # pkg = ".",
    document = TRUE,
    # build_args = NULL,
    # manual = FALSE,
    # cran = TRUE,
    # remote = FALSE,
    # force_suggests = FALSE,
    run_dont_test = TRUE,
    args = "--timings",
    # env_vars = c(NOT_CRAN = "true"),
    # quiet = FALSE,
    # check_dir = NULL,
    # cleanup = deprecated(),
    vignettes = FALSE,
    error_on = "warning"  # c("never", "error", "warning", "note")
)
