* Despite repeated efforts we did not manage to replicate the clang-UBSAN errors detected on CRAN.
* We added additional checks in the compiled code which should ensure that the pointed out undefined behaviour error is avoided, but due to not being able to replicate the error, we could not test if this fixes the error.
