test_splitIndicies <- function()
{
    .splitIndices <- BiocParallel:::.splitIndices

    checkIdentical(list(), .splitIndices(0, 0))
    checkIdentical(list(), .splitIndices(0, 1))
    checkIdentical(list(), .splitIndices(0, 2))

    checkIdentical(list(1:4), .splitIndices(4, 0))
    checkIdentical(list(1:4), .splitIndices(4, 1))
    checkIdentical(list(1:2, 3:4), .splitIndices(4, 2))
    checkIdentical(as.list(1:4), .splitIndices(4, 4))
    checkIdentical(as.list(1:4), .splitIndices(4, 8))

    checkIdentical(list(1:4, 5:7), .splitIndices(7, 2))
}

test_splitX <- function()
{
    .splitX <- BiocParallel:::.splitX

    checkIdentical(list(), .splitX(character(), 0, 0))
    checkIdentical(list(), .splitX(character(), 1, 0))
    checkIdentical(list(), .splitX(character(), 0, 1))
    checkIdentical(list(), .splitX(character(), 1, 1))

    X <- LETTERS[1:4]
    checkIdentical(list(X),              .splitX(X, 0, 0))
    checkIdentical(list(X),              .splitX(X, 1, 0))
    checkIdentical(list(X[1:2], X[3:4]), .splitX(X, 2, 0))
    checkIdentical(as.list(X),           .splitX(X, 4, 0))
    checkIdentical(as.list(X),           .splitX(X, 8, 0))

    checkIdentical(list(X[1:2], X[3:4]), .splitX(X, 2, 0))
    checkIdentical(list(X), .splitX(X, 2, 1))
    checkIdentical(list(X[1:2], X[3:4]),     .splitX(X, 2, 2))
    checkIdentical(list(X[1], X[2:3], X[4]), .splitX(X, 2, 3))
    checkIdentical(as.list(X),               .splitX(X, 2, 4))
}

test_splitX_with_redo <- function(){
    .split_X_redo <- BiocParallel:::.split_X_redo
    .bploop_rng_iter <- BiocParallel:::.bploop_rng_iter

    N <- 6
    X <- 1:N


    task_size <- 1L
    redo_index <- rep(TRUE, N)
    n <- .split_X_redo(X[redo_index], redo_index, task_size)
    splittedX <- .split_X_redo(X[redo_index], redo_index, task_size, n)
    checkIdentical(splittedX, as.list(X))

    task_size <- 4L
    redo_index <- rep(TRUE, N)
    n <- .split_X_redo(X[redo_index], redo_index, task_size)
    splittedX <- .split_X_redo(X[redo_index], redo_index, task_size, n)
    checkIdentical(splittedX, list(1:4, 5:6))

    task_size <- 1L
    redo_index <- c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE)
    n <- .split_X_redo(X[redo_index], redo_index, task_size)
    splittedX <- .split_X_redo(X[redo_index], redo_index, task_size, n)
    checkIdentical(
        splittedX,
        list(1L, .bploop_rng_iter(1L),
             3L, 4L, 5L)
    )

    task_size <- 2L
    redo_index <- c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE)
    n <- .split_X_redo(X[redo_index], redo_index, task_size)
    splittedX <- .split_X_redo(X[redo_index], redo_index, task_size, n)
    checkIdentical(
        splittedX,
        list(1L, .bploop_rng_iter(1L),
             3L:4L, 5L)
    )

    task_size <- 10L
    redo_index <- c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE)
    n <- .split_X_redo(X[redo_index], redo_index, task_size)
    splittedX <- .split_X_redo(X[redo_index], redo_index, task_size, n)
    checkIdentical(
        splittedX,
        list(1L, .bploop_rng_iter(1L),
             3L:5L)
    )

    task_size <- 1L
    redo_index <- c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)
    n <- .split_X_redo(X[redo_index], redo_index, task_size)
    splittedX <- .split_X_redo(X[redo_index], redo_index, task_size, n)
    checkIdentical(
        splittedX,
        list(.bploop_rng_iter(2L),
             3L)
    )

    task_size <- 1L
    redo_index <- c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE)
    n <- .split_X_redo(X[redo_index], redo_index, task_size)
    splittedX <- .split_X_redo(X[redo_index], redo_index, task_size, n)
    checkIdentical(
        splittedX,
        list(.bploop_rng_iter(2L),
             3L, 4L, 5L, 6L)
    )

    task_size <- 3L
    redo_index <- c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE)
    n <- .split_X_redo(X[redo_index], redo_index, task_size)
    splittedX <- .split_X_redo(X[redo_index], redo_index, task_size, n)
    checkIdentical(
        splittedX,
        list(.bploop_rng_iter(2L),
             3:5, 6L)
    )
}

test_redo_index <- function() {
    .redo_index <- BiocParallel:::.redo_index
    err <- BiocParallel:::.error("")
    checkIdentical(logical(), .redo_index(list(), list()))
    checkIdentical(TRUE, .redo_index(list(1), list(err)))
    checkIdentical(c(FALSE, TRUE), .redo_index(list(1, "2"), list(1, err)))
    ## all need recalculating
    checkIdentical(c(TRUE, TRUE), .redo_index(list("1", "2"), list(err, err)))
    ## X can be a vector
    checkIdentical(c(FALSE, TRUE), .redo_index(1:2, list(1, err)))
    ## lengths differ
    checkException(.redo_index(list(1, 2), list(err)), silent=TRUE)
    ## no previous error
    checkException(.redo_index(list(1, 2), list(1, 2)), silent=TRUE)
}

test_rename <- function() {
    .rename <- BiocParallel:::.rename

    X <- list()
    Y <- character()
    Z <- list(X)
    W <- list(Y)
    checkIdentical(X, .rename(list(), X))
    checkIdentical(X, .rename(list(), Y))
    checkIdentical(X, .rename(list(), Z))
    checkIdentical(X, .rename(list(), W))

    names(X) <- names(Y) <- character()
    Z <- list(X)
    W <- list(Y)
    checkIdentical(X, .rename(list(), X))
    checkIdentical(X, .rename(list(), Y))
    checkIdentical(list(), .rename(list(), Z))
    checkIdentical(list(), .rename(list(), W))

    Z <- list(x = X)
    W <- list(x = Y)
    checkIdentical(list(x = 1), .rename(list(1), Z))
    checkIdentical(list(x = 1), .rename(list(1), W))

    X <- list(a = 1:2)
    exp0 <- vector("list", length(X))
    checkIdentical(setNames(exp0, names(X)), .rename(exp0, X))

    X <- list(c(a = 1))
    exp0 <- vector("list", length(X))
    checkIdentical(exp0, .rename(exp0, X))

    Y <- c(x = "a")
    checkIdentical(Y, .rename(Y, Y))

    X <- list(a = 1:2, b = 3:4)
    exp0 <- vector("list", length(X))
    exp <- setNames(exp0, names(X))
    checkIdentical(exp, .rename(exp0, X))

    X <- list(c(a = 1))
    exp0 <- vector("list", length(X))
    checkIdentical(exp0, .rename(exp0, X))

    X <- list(A = c(a = 1, b=2), B = c(c = 1, d = 2))
    exp0 <- vector("list", length(X))
    exp <- setNames(exp0, names(X))
    checkIdentical(exp, .rename(exp0, X))
}

test_mrename <- function() {
    .mrename <- BiocParallel:::.mrename

    ## list() when `mapply()` invoked with no arguments, `mapply(identity)`
    X <- list()
    checkIdentical(X, .mrename(list(), X))

    ## list(X) when `mapply()` invoked with one argument, `mapply(identity, X)`
    XX <- list(X)
    checkIdentical(list(), .mrename(list(), XX))

    ## `mapply(identity, character())` returns a _named_ list()
    X <- character()
    XX <- list(X)
    checkIdentical(setNames(list(), character()), .mrename(list(), XX))

    ## named arguments to mapply() are _not_ names of return value...
    X = list()
    XX <- list(x = X)
    checkIdentical(list(), .mrename(list(), XX))

    ## ...except if the argument is a character()
    X <- character()
    XX <- list(x = X)
    checkIdentical(setNames(list(), character()), .mrename(list(), XX))

    ## with multiple arguments, names are from the first argument
    XX = list(c(a = 1, b = 2, c = 3), c(d = 4, e = 5, f = 6))
    checkIdentical(
        setNames(list(1, 2, 3), letters[1:3]),
        .mrename(list(1, 2, 3), XX)
    )

    ## ...independent of names on the arguments
    XX = list(A = c(a = 1, b = 2, c = 3), B = c(d = 4, e = 5, f = 6))
    checkIdentical(
        setNames(list(1, 2, 3), letters[1:3]),
        .mrename(list(1, 2, 3), XX)
    )

    ## when the first argument is an unnamed character vector, names are values
    XX = list(A = letters[1:3], B = 1:3)
    checkIdentical(
        setNames(list(1, 2, 3), letters[1:3]),
        .mrename(list(1, 2, 3), XX)
    )

    ## ...except if there are names on the first vector...
    XX = list(A = setNames(letters[1:3], LETTERS[1:3]), B = 1:3)
    checkIdentical(
        setNames(list(1, 2, 3), LETTERS[1:3]),
        .mrename(list(1, 2, 3), XX)
    )
}
