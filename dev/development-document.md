# Development Document

## Build for CRAN

## Build the package

```sh
# Create documentation (optional)
./dev/documentation-build.sh
# Build the package
./dev/run-cran-build.sh
# Test the package for CRAN
./dev/run-cran-check.sh
```

## Scripts in dev folder

- [documentation-build.R](./documentation-build.R): Script to build the documentation.
- [package-setup.R](./package-setup.R): Script to setup a new package.
- [run-check.R](./run-check.R): Script to check the package.
- [run-cran-build.sh](./run-cran-build.sh): Script to build the package.
- [run-cran-check.sh](./run-cran-check.sh): Script to test the package for CRAN.
- [run-test.R](./run-test.R): Script to run the unit tests.