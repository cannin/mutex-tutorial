#' Extracts a three column edge table from correlation/adjacency matrix with top interactions
#'
#' @param input a named symmetric square matrix with pair-wise correlation values
#' both colnames and rownames should be the same, a vector of strings of node names
#'
#' @param n a numeric value of how many top interactions to choose
#' 
#' @param phenoName a vector of strings of the phenotypic nodes
#' 
#' @param actName a vector of strings of the activity/drug coupling nodes
#'
#' @return extractEdge returns a three column matrix: participant A, participant B, partial correlation
#'
#' @note phenoName and actName are optional arguments but 
#' make sure the phenoName and actName are subsets of input colnames/rownames
#' 
#' @examples
#' extractEdge(inputMatrix, 150, c("G1arrest","G2arrest","cellviab"), c("aMEK","aAKT","aHDAC"))
#'
#' @author Judy Shen (c_shen@g.harvard.edu)
extractEdge <- function(input, n=NULL, phenoName, actName, fdrMat=NULL) {
  edges_list <- data.frame(PARTICIPANT_A=character(length=2*n),
                           PARTICIPANT_B=character(length=2*n),
                           PARTIAL_CORRELATION=numeric(length=2*n),
                           FDR_VAL=numeric(length=2*n),
                           stringsAsFactors = FALSE)
  # ignore self-interaction, drug-drug interaction, phenotype-phenotype interaction
  # and direct drug-phenotype interaction
  diag(input) <- 0
  
  if (!missing(phenoName) & !missing(actName)) {
    input[c(phenoName, actName), c(phenoName,actName)] <- 0
  } else if (!missing(phenoName) & missing(actName)) {
    input[phenoName, phenoName] <- 0
  } else if (missing(phenoName) & !missing(actName)) {
    input[actName, actName] <- 0
  }
  
  top_edges <- order(abs(input), decreasing = TRUE)[1:(2*n)]
  
  for (i in 1:(2*n)) {
    idx <- arrayInd(top_edges[i],dim(input))
    
    if(!is.null(fdrMat)) {
      edges_list[i,] <- c(colnames(input)[idx[1,2]], rownames(input)[idx[1,1]], input[idx[1,1],idx[1,2]], fdrMat[idx[1,1],idx[1,2]])
    } else {
      edges_list[i,] <- c(colnames(input)[idx[1,2]], rownames(input)[idx[1,1]], input[idx[1,1],idx[1,2]])  
    }
  }
  
  return(edges_list[1:nrow(edges_list)%%2==0,])
}
