exon_overlap <- function(gene_A_exons, gene_B_exons){
  
  if(dim(gene_A_exons)[1] == 0 | dim(gene_B_exons)[1] == 0){
    return (FALSE)
  }

  for(row_exonA in 1:nrow(gene_A_exons)){
    for(row_exonB in 1:nrow(gene_B_exons)){

      if (any(gene_A_exons[,1] <= 0 | gene_A_exons[,2] <= 0) | 
        any(gene_B_exons[,1] <= 0 | gene_B_exons[,2] <= 0)) {
        return(TRUE)  
      }

      from_x <- gene_A_exons[row_exonA,1]
      to_x <- gene_A_exons[row_exonA,2]
      if (to_x > from_x) {
        to_x <- to_x - 1
      }
      x = seq(from = from_x, to = to_x, by = 1)

      from_y <- gene_B_exons[row_exonB,1]
      to_y <- gene_B_exons[row_exonB,2]
      if (to_y > from_y) {
        to_y <- to_y - 1
      }
      y = seq(from = from_y, to = to_y, by = 1)

      if(length(intersect(x,y))!=0){
        return (TRUE)
      }
    }
  }
  return (FALSE)
}