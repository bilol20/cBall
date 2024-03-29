# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

CppS <- function(n, Kz, Wz, Kyz, Wyz, L) {
    .Call(`_cBall_CppS`, n, Kz, Wz, Kyz, Wyz, L)
}

rowcolsampler <- function(A, s) {
    .Call(`_cBall_rowcolsampler`, A, s)
}

resample <- function(n, Kz, Wz, Kyz, Wyz, L, Pi) {
    .Call(`_cBall_resample`, n, Kz, Wz, Kyz, Wyz, L, Pi)
}

kernel <- function(x, h) {
    .Call(`_cBall_kernel`, x, h)
}

fun2 <- function(x) {
    .Call(`_cBall_fun2`, x)
}

eu <- function(x, y) {
    .Call(`_cBall_eu`, x, y)
}

dist_cpp <- function(X) {
    .Call(`_cBall_dist_cpp`, X)
}

