
names(mtcars)
two_part = function(formula = mpg ~ cyl + hp + wt | cyl + am, data = mtcars) { 

  ## extracted from pscl::zero_inflation, identifies "|" in formulas. 
    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)#;  return(mf)
    
    m <- match(c("formula", "data", "subset", "na.action", "weights", 
        "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    
    if ( length(formula[[3]]) > 1 && 
         identical(formula[[3]][[1]], as.name("|")) ) {
        ff <- formula
        formula[[3]][1] <- call("+")
        mf$formula <- formula
        ffc <- . ~ .
        ffz <- ~.
        ffc[[2]] <- ff[[2]]
        ffc[[3]] <- ff[[3]][[2]]
        ffz[[3]] <- ff[[3]][[3]]
        ffz[[2]] <- NULL
    } else {
        ffz <- ffc <- ff <- formula
        ffz[[2]] <- NULL
    }
    if (inherits(try(terms(ffz), silent = TRUE), "try-error")) {
        ffz <- eval(parse(text = sprintf(paste("%s -", deparse(ffc[[2]])), 
            deparse(ffz))))
    }
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    mtX <- terms(ffc, data = data)
    X <- model.matrix(mtX, mf)
    mtZ <- terms(ffz, data = data)
    mtZ <- terms(update(mtZ, ~.), data = data)
    Z <- model.matrix(mtZ, mf)
    Y <- model.response(mf, "numeric")
  ## end code taken from pscl::zero_inflation
     
    list( Y = Y, X = X, Z = Z)
}

test0 = two_part(formula = mpg ~ cyl + hp + wt | cyl + am, data = mtcars)
test1 = two_part(formula = mpg ~ cyl + hp + wt, data = mtcars)
stopifnot( all.equal( test1$X, test1$Z) )
