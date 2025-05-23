test_rng_fun_advances_generator <- function()
{
    .rng_get_generator <- BiocParallel:::.rng_get_generator
    .rng_reset_generator <- BiocParallel:::.rng_reset_generator
    .rng_job_fun_factory <- BiocParallel:::.rng_job_fun_factory
    .RNGstream <- BiocParallel:::.RNGstream

    state <- .rng_get_generator()
    on.exit(.rng_reset_generator(state$kind, state$seed))

    SEED <- .RNGstream(bpstart(SerialParam()))
    checkIdentical(
        ## independently invoked with same seed --> same result
        .rng_job_fun_factory(function(i) rnorm(i), SEED)(2),
        .rng_job_fun_factory(function(i) rnorm(i), SEED)(2)
    )

    SEED <- .RNGstream(bpstart(SerialParam()))
    FUN <- .rng_job_fun_factory(function(i) rnorm(i), SEED)
    target <- FUN(2) # two numbers from same stream

    FUN <- .rng_job_fun_factory(function(i) rnorm(i), SEED)
    obs <- c(FUN(1), FUN(1)) # two numbers from separate streams
    checkIdentical(obs[[1]], target[[1]])
    checkTrue(obs[[2]] != target[[2]])

    checkTrue(identical(state, .rng_get_generator()))
}

test_rng_lapply <- function()
{
    .rng_get_generator <- BiocParallel:::.rng_get_generator
    .rng_reset_generator <- BiocParallel:::.rng_reset_generator

    .RNGstream <- BiocParallel:::.RNGstream
    .rng_lapply <- BiocParallel:::.rng_lapply
    .rng_next_substream <- BiocParallel:::.rng_next_substream

    state <- .rng_get_generator()
    on.exit(.rng_reset_generator(state$kind, state$seed))

    SEED <- .RNGstream(bpstart(SerialParam()))
    checkIdentical(
        ## same sequence of random number streams
        .rng_lapply(1:2, function(i) rnorm(1), BPRNGSEED = SEED),
        .rng_lapply(1:2, function(i) rnorm(1), BPRNGSEED = SEED)
    )

    SEED1 <- .RNGstream(bpstart(SerialParam()))
    SEED2 <- .rng_next_substream(SEED1)
    target <- .rng_lapply(1:2, function(i) rnorm(2), BPRNGSEED = SEED1)
    obs <- c(
        .rng_lapply(1, function(i) rnorm(2), BPRNGSEED = SEED1),
        .rng_lapply(1, function(i) rnorm(2), BPRNGSEED = SEED2)
    )
    checkIdentical(target, obs)

    checkTrue(identical(state, .rng_get_generator()))
}

test_rng_bplapply <- function()
{
    .rng_get_generator <- BiocParallel:::.rng_get_generator
    .rng_reset_generator <- BiocParallel:::.rng_reset_generator

    state <- .rng_get_generator()
    on.exit(.rng_reset_generator(state$kind, state$seed))

    p1 <- SerialParam(RNGseed = 123)
    p2 <- SnowParam(3, RNGseed = 123)
    p3 <- SnowParam(5, RNGseed = 123)
    FUN <- function(i) rnorm(2)

    ## SerialParam / SnowParam same results
    target <- bplapply(1:11, FUN, BPPARAM = p1)
    checkIdentical(bplapply(1:11, FUN, BPPARAM = p2), target)

    ## SerialParam / SnowParam same results, different number of workers
    checkIdentical(bplapply(1:11, FUN, BPPARAM = p3), target)

    if (identical(.Platform$OS.type, "unix")) {
        ## SerialParam / TransientMulticoreParam same results
        p4a <- MulticoreParam(5, RNGseed = 123)
        checkIdentical(bplapply(1:11, FUN, BPPARAM = p4a), target)
        ## SerialParam / MulticoreParam same results
        p4b <- bpstart(MulticoreParam(5, RNGseed = 123))
        checkIdentical(bplapply(1:11, FUN, BPPARAM = p4b), target)
        bpstop(p4b)
    }

    ## single worker coerced to SerialParam
    p5 <- SnowParam(1, RNGseed = 123)
    checkIdentical(bplapply(1:11, FUN, BPPARAM = p5), target, "p5")
    checkIdentical(state$kind, .rng_get_generator()$kind)
}

test_rng_bpiterate <- function()
{
    .rng_get_generator <- BiocParallel:::.rng_get_generator
    .rng_reset_generator <- BiocParallel:::.rng_reset_generator

    state <- .rng_get_generator()
    on.exit(.rng_reset_generator(state$kind, state$seed))

    FUN <- function(i) rnorm(2)

    ITER_factory <- function() {
        x <- 1:11
        i <- 0L
        function() {
            i <<- i + 1L
            if (i > length(x))
                return(NULL)
            x[[i]]
        }
    }

    p1 <- SerialParam(RNGseed = 123)
    p2 <- SnowParam(3, RNGseed = 123)
    p3 <- SnowParam(5, RNGseed = 123)

    target <- bplapply(1:11, FUN, BPPARAM = p1)
    checkIdentical(target, bpiterate(ITER_factory(), FUN, BPPARAM = p1), "p1")
    checkIdentical(target, bpiterate(ITER_factory(), FUN, BPPARAM = p2), "p2")
    checkIdentical(target, bpiterate(ITER_factory(), FUN, BPPARAM = p3), "p3")

    if (identical(.Platform$OS.type, "unix")) {
        ## SerialParam / TransientMulticoreParam same results
        p4a <- MulticoreParam(5, RNGseed = 123)
        checkIdentical(
            target, bpiterate(ITER_factory(), FUN, BPPARAM = p4a),
            "p4a"
        )
        ## SerialParam / MulticoreParam same results
        p4b <- bpstart(MulticoreParam(5, RNGseed = 123))
        checkIdentical(
            target, bpiterate(ITER_factory(), FUN, BPPARAM = p4b),
            "p4b"
        )
        bpstop(p4b)
    }
    checkIdentical(state$kind, .rng_get_generator()$kind)
}

test_rng_bpstart <- function()
{
    .rng_get_generator <- BiocParallel:::.rng_get_generator
    .rng_reset_generator <- BiocParallel:::.rng_reset_generator
    state <- .rng_get_generator()

    FUN <- function(i) rnorm(2)

    ITER_factory <- function() {
        x <- 1:11
        i <- 0L
        function() {
            i <<- i + 1L
            if (i > length(x))
                return(NULL)
            x[[i]]
        }
    }

    ## bplapply
    p0 <- bpstart(SerialParam()) # random seed
    result1 <- unlist(bplapply(1:11, FUN, BPPARAM = p0))
    result2 <- unlist(bplapply(1:11, FUN, BPPARAM = p0))
    checkTrue(!any(result1 %in% result2))

    bpstart(bpstop(p0)) # different random seed
    result3 <- unlist(bplapply(1:11, FUN, BPPARAM = p0))
    checkTrue(!any(result3 %in% result1))

    p0 <- bpstart(SerialParam(RNGseed = 123)) # set seed
    result1 <- unlist(bplapply(1:11, FUN, BPPARAM = p0))
    result2 <- unlist(bplapply(1:11, FUN, BPPARAM = p0))
    checkTrue(!any(result1 %in% result2))

    bpstart(bpstop(p0)) # reset seed, same stream
    result3 <- unlist(bplapply(1:11, FUN, BPPARAM = p0))
    result4 <- unlist(bplapply(1:11, FUN, BPPARAM = p0))
    checkIdentical(result3, result1)
    checkIdentical(result4, result2)

    ## bpiterate
    p0 <- bpstart(SerialParam()) # random seed
    result1 <- unlist(bpiterate(ITER_factory(), FUN, BPPARAM = p0))
    result2 <- unlist(bpiterate(ITER_factory(), FUN, BPPARAM = p0))
    checkTrue(!any(result1 %in% result2))

    bpstart(bpstop(p0)) # different random seed
    result3 <- unlist(bpiterate(ITER_factory(), FUN, BPPARAM = p0))
    checkTrue(!any(result3 %in% result1))

    p0 <- bpstart(SerialParam(RNGseed = 123)) # set seed
    result1 <- unlist(bpiterate(ITER_factory(), FUN, BPPARAM = p0))
    result2 <- unlist(bpiterate(ITER_factory(), FUN, BPPARAM = p0))
    checkTrue(!any(result1 %in% result2))

    bpstart(bpstop(p0)) # reset seed, same stream
    result3 <- unlist(bpiterate(ITER_factory(), FUN, BPPARAM = p0))
    result4 <- unlist(bpiterate(ITER_factory(), FUN, BPPARAM = p0))
    checkIdentical(result3, result1)
    checkIdentical(result4, result2)

    checkIdentical(state$kind, .rng_get_generator()$kind)
}

.test_rng_bpstart_does_not_iterate_rng_seed <- function(param) {
    .rng_get_generator <- BiocParallel:::.rng_get_generator

    state <- .rng_get_generator()
    set.seed(123L)
    target <- runif(1L)

    ## bpstart() should not increment the random number seed
    set.seed(123L)
    bpstart(param)
    checkIdentical(target, runif(1L))
    bpstop(param)

    ## bplapply does not increment stream
    set.seed(123)
    result <- bplapply(1:3, runif, BPPARAM = param)
    checkIdentical(target, runif(1L))

    ## bplapply uses internal stream
    set.seed(123)
    result <- bplapply(1:3, runif, BPPARAM = param)
    checkTrue(!identical(result, bplapply(1:3, runif, BPPARAM = param)))
    checkIdentical(target, runif(1L))
    target1 <- lapply(1:3, runif)
    checkTrue(!identical(result, target1))

    checkIdentical(state$kind, .rng_get_generator()$kind)
}

test_rng_bpstart_does_not_iterate_rng_seed <- function() {
    .rng_get_generator <- BiocParallel:::.rng_get_generator
    .rng_reset_generator <- BiocParallel:::.rng_reset_generator
    TEST_FUN <- .test_rng_bpstart_does_not_iterate_rng_seed

    state <- .rng_get_generator()
    on.exit(.rng_reset_generator(state$kind, state$seed))

    TEST_FUN(SerialParam())
    TEST_FUN(SnowParam(2))
    if (identical(.Platform$OS.type, "unix"))
        TEST_FUN(MulticoreParam(2))
}

.test_rng_global_and_RNGseed_indepenent <- function(param_fun) {
    set.seed(123)
    target <- bplapply(1:3, runif, BPPARAM = param_fun())
    current <- bplapply(1:3, runif, BPPARAM = param_fun(RNGseed = 123))
    checkTrue(!identical(target, current))
}

test_rng_global_and_RNGseed_independent <- function() {
    .rng_get_generator <- BiocParallel:::.rng_get_generator
    .rng_reset_generator <- BiocParallel:::.rng_reset_generator
    TEST_FUN <- .test_rng_global_and_RNGseed_indepenent

    state <- .rng_get_generator()
    on.exit(.rng_reset_generator(state$kind, state$seed))

    TEST_FUN(SerialParam)
    TEST_FUN(SnowParam)
    if (identical(.Platform$OS.type, "unix"))
        TEST_FUN(MulticoreParam)
}

.test_rng_bpredo_impl <- function(param) {
    FUN <- function(i) rnorm(1)
    target <- unlist(bplapply(1:11, FUN, BPPARAM = param))

    FUN0 <- function(i) {
        if (identical(i, 7L)) {
            stop("i == 7")
        } else rnorm(1)
    }
    result <- bptry(bplapply(1:11, FUN0, BPPARAM = param))
    checkIdentical(unlist(result[-7]), target[-7])
    checkTrue(inherits(result[[7]], "remote_error"))

    FUN1 <- function(i) {
        if (identical(i, 7L)) {
            ## the random number stream should be in the same state as the
            ## first time through the loop, and rnorm(1) should return
            ## same result as FUN
            rnorm(1)
        } else {
            ## if this branch is used, then we are incorrectly updating
            ## already calculated elements -- '0' in the output would
            ## indicate this error
            0
        }
    }
    result <- unlist(bplapply(1:11, FUN1, BPREDO = result, BPPARAM = param))
    checkIdentical(result, target)


    bpstart(param)
    target1 <- unlist(bplapply(1:11, FUN, BPPARAM = param))
    target2 <- unlist(bplapply(1:11, FUN, BPPARAM = param))
    target3 <- unlist(bplapply(1:11, FUN, BPPARAM = param))
    bpstop(param)

    bpstart(param)
    result1 <- bptry(bplapply(1:11, FUN0, BPPARAM = param))
    result1_redo1 <- unlist(bplapply(1:11, FUN1, BPREDO = result1, BPPARAM = param))
    result2 <- unlist(bplapply(1:11, FUN, BPPARAM = param))
    result1_redo2 <- unlist(bplapply(1:11, FUN1, BPREDO = result1, BPPARAM = param))
    result3 <- unlist(bplapply(1:11, FUN, BPPARAM = param))
    checkIdentical(target1, result1_redo1)
    checkIdentical(target1, result1_redo2)
    checkIdentical(target2, result2)
    checkIdentical(target3, result3)
}

test_rng_bpredo <- function()
{
    .rng_get_generator <- BiocParallel:::.rng_get_generator
    .rng_reset_generator <- BiocParallel:::.rng_reset_generator

    state <- .rng_get_generator()
    on.exit(.rng_reset_generator(state$kind, state$seed))

    param <- SerialParam(RNGseed = 123, stop.on.error = FALSE)
    .test_rng_bpredo_impl(param)

    if (identical(.Platform$OS.type, "unix")) {
        param <- MulticoreParam(3, RNGseed = 123, stop.on.error = FALSE)
        .test_rng_bpredo_impl(param)
    }
}

test_rng_fallback_SerialParam <- function()
{
    .rng_get_generator <- BiocParallel:::.rng_get_generator
    .rng_reset_generator <- BiocParallel:::.rng_reset_generator

    state <- .rng_get_generator()
    on.exit(.rng_reset_generator(state$kind, state$seed))

    FUN <- function(i) rnorm(1)
    param <- SerialParam(RNGseed = 123, stop.on.error = FALSE)
    target <- unlist(bplapply(1:2, FUN, BPPARAM = param))

    param2 <- SnowParam(RNGseed = 123, workers = 1L)
    checkIdentical(unlist(bplapply(1:2, FUN, BPPARAM = param2)), target)

    bpstart(param2)
    checkIdentical(bplapply(1, FUN, BPPARAM = param2)[[1]], target[1])
    checkIdentical(bplapply(1, FUN, BPPARAM = param2)[[1]], target[2])
    bpstop(param2)
}
