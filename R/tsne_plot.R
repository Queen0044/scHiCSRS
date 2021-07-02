#' tsne_plot
#'
#' @param data The observed matrix, with each column being the uppertriangular of a single cell HiC matrix.
#' @param cell_type A vector that indicates cell type.
#' @param dims Integer. Output dimentionality. Default=2.
#' @param perplexity Numeric; Perplexity parameter (should not be bigger than 3 * perplexity < nrow(X) - 1, see details for interpretation).
#' @param check_duplicates Logical; Checks whether duplicates are present. It is best to make sure there are no duplicates present and set this option to FALSE, especially for large datasets (default: TRUE).
#' @param seed Random seed.
#' @param kmeans Logical. Whether apply kmeans clustering on tsne data or not.
#' @param ncenters Number of clusters in kmeans clustering.
#'
#' @return
#' A stne visualization plot.
#' @export
#'
#' @import Rtsne
#'
#' @import ggpubr
#'
#' @import ggplot2
#'
#' @examples
#' tsne_plot(simudat, cell_type=c(rep("L4", 131), rep("L5",180)), seed=1250)
tsne_plot <- function(data, cell_type, dims = 2, perplexity=10, check_duplicates = FALSE, seed=1234,
                      kmeans=TRUE, ncenters) {

  mydata <- t(scale(data)) # standardize variables
  #rownames(mydata) = c(paste("A",1:131,sep="_"), paste("B",1:180,sep="_"))

  set.seed(seed)
  tsne_dat <- Rtsne(mydata, perplexity=perplexity, check_duplicates = check_duplicates)

  data=tsne_dat
  if(kmeans){
    my.xy = data.frame(x=data$Y[,1], y=data$Y[,2])

    res.km <- kmeans(my.xy, centers = ncenters, nstart = 25)

    my.xy$cluster <- factor(res.km$cluster)
    my.xy$type <- cell_type

    ggscatter(
      my.xy, x = "x", y = "y", xlab = "",
      ylab = "", color = "cluster", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
      shape = "type", size = 2,  legend = "right", ggtheme = theme_bw()) +
      scale_shape_manual(values=c(2, 19, 17, 10)[1:ncenters])
  }else{
    my.xy = data.frame(x=data$Y[,1], y=data$Y[,2], group=cell_type)
    ggplot(my.xy) + geom_point(aes(x=x, y=y, color=group))
  }
}






