# bios_1000pts

Runs the final step (step 6. Final centennial run) of the BIOS 1000pts
configuration initialised from restart files.

## How to generate a code coverage report:

1. Add `-prof-gen=srcpos -prof-dir=<path/to/profiling/output/dir>` to Intel
   compiler flags.
2. Compile and run CABLE.
3. [Gadi] Load the Intel compiler module: `module load intel-compiler`.
4. Change directory into the profiling output directory and run `profmerge -a`.
5. Generate HTML coverage report: `codecov -prj myProject -spi pgopti.spi -dpi
   pgopti.dpi`. This creates a single HTML file (CODE_COVERAGE.HTML) and a
   sub-directory (CodeCoverage) in the current directory. Open the file in a
   web browser to view the reports.

See [Code Coverage Tool][code-coverage-tool] for more information.

[code-coverage-tool]: https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2023-0/code-coverage-tool.html
