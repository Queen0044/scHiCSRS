#' Clean the data
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
#' clean.data(X_count)
clean.data <- function(x) {

  if (!(grepl("matrix", class(x), ignore.case = TRUE))) {

    x <- Matrix::Matrix(as.matrix(x))

    message("Converting x to matrix.")

    if (!is.numeric(x)) {

      warning("Make sure x is numeric.")

    }

  }

  np <- dim(x)

  size <- as.numeric(np[1])*as.numeric(np[2])

  if(size > 2^31-1){

    inds <- split(1:np[2], ceiling(1:np[2]/1000))

    for(i in 1:length(inds)){

      x[, inds[[i]]][x[, inds[[i]]] < 0.001] <- 0

    }

  } else {

    x[x < 0.001] <- 0

  }

  if (is.null(np) | (np[2] <= 1))

    stop("x should be a matrix with 2 or more columns")

  if (min(Matrix::colSums(x)) == 0) {

    nzerocells <- sum(Matrix::colSums(x) == 0)

    x <- x[, Matrix::colSums(x) != 0]

    message("Removing ", nzerocells, " cell(s) with zero expression.")

  }

  if (is.null(rownames(x))) {

    rownames(x) <- 1:np[1]

  }

  x

}
