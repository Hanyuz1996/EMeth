pkgname <- "EMeth"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('EMeth')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("EMeth-package")
### * EMeth-package

flush(stderr()); flush(stdout())

### Name: EMeth-package
### Title: Cell type decomposition based on methylation data
### Aliases: EMeth-package EMeth

### ** Examples

# Examples can be found in the "example" subdirectory.



cleanEx()
nameEx("cv.emeth")
### * cv.emeth

flush(stderr()); flush(stdout())

### Name: cv.emeth
### Title: cv.emeth: cross validation for emeth.
### Aliases: cv.emeth

### ** Examples

## See examples folder.



cleanEx()
nameEx("emeth")
### * emeth

flush(stderr()); flush(stdout())

### Name: emeth
### Title: emeth: cell type decomposition from DNA methylation data based
###   on EM-type algorithm and penalized likelihood maximization.
### Aliases: emeth

### ** Examples

## See examples folder



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
