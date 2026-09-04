# Created by Peter Crisp
# August 2014
# pedrocrisp at gmail.com

# For use in differential gene expression analysis.

##### MAplotGeneSet function #####
  
  #for debugging
  # input_list = "xrn2xrn3.vs.control"
  # input_list_name = "xrn2xrn3"
  # gene_set_list = "alx8.vs.control"
  # gene_set_list_name = "alx8"
  # marray_lmobject = fit2

MAplotGeneSetLimma <-  function(marray_lmobject,
                                results_matrix,
                                input_list,
                                input_list_name,
                                gene_set_list,
                                gene_set_list_name){
  limma::plotMA(marray_lmobject,
                coef   = input_list,
                status = results_matrix[, gene_set_list],
                cex    = 0.5, values = c(1, -1), col = c("red", "blue"),
                main   = paste0(input_list), legend = FALSE)
  abline(h = 0, col = "darkgrey")
  legend("bottomright",
         legend = (c(paste0("up in ", gene_set_list_name),
                     paste0("down in ", gene_set_list_name))),
         bty =  "n", pch = c(16, 16),
         col = c("red", "blue"), text.col = c("red", "blue"))
}

ma_plot_geneset_green_limma <- function(marray_lmobject,
                                        results_matrix,
                                        input_list,
                                        input_list_name,
                                        gene_set_list,
                                        gene_set_list_name){

  limma::plotMA(marray_lmobject,
                coef   = input_list,
                status = results_matrix[, gene_set_list],
                cex    = 0.5, values = c(1, -1),
                col    = c("red", "green"),
                main   = paste0(input_list), legend = FALSE)
  abline(h = 0, col = "darkgrey")
  legend("bottomright",
         legend =  (c(paste0("up in ",   gene_set_list_name),
                      paste0("down in ", gene_set_list_name))),
         bty    =  "n", pch = c(16, 16),
         col    = c("red", "green"), text.col = c("red", "green"))
}

ma_plot_geneset_2limma <- function(marray_lmobject,
                                   geneset_results_matrix,
                                   de_call_names,
                                   input_list, input_list_name,
                                   gene_set_list_name) {
  limma::plotMA(marray_lmobject,
                coef   = input_list,
                status = geneset_results_matrix,
                cex    = 0.35, values = de_call_names,
                col    = c("#00FF00", "#FF0000"),
                main   = paste0(input_list_name), legend = FALSE)
  abline(h = 0, col = "darkgrey")
  legend("bottomright",
         legend   = (c(paste0("up in ", gene_set_list_name),
                       paste0("down in ", gene_set_list_name))),
         cex      = 0.5, bty = "n", pch = c(16, 16),
         col      = c("#00FF00", "#FF0000"),
         text.col = c("black",   "black"))
}

ma_plot_geneset_2r_limma <-  function(marray_lmobject,
                                      geneset_results_matrix,
                                      de_call_names,
                                      input_list, input_list_name,
                                      gene_set_list_name) {
  limma::plotMA(marray_lmobject,
                coef   = input_list,
                cex    = 0.35, col = c("red", "blue"),
                status = geneset_results_matrix, values = de_call_names,
                main   = paste0(input_list_name), legend = FALSE)

  abline(h = 0, col = "darkgrey")
  legend("bottomright",
         cex      = 0.5, col = c("red", "blue"),
         text.col = c("black", "black"),
         legend   =  (c(paste0("up in ",   gene_set_list_name),
                        paste0("down in ", gene_set_list_name))),
         bty      = "n", pch = c(16, 16))
}

MAplotGeneSet2gLimma <- function(marray_lmobject,
                                 geneset_results_matrix,
                                 de_call_names,
                                 input_list, input_list_name,
                                 gene_set_list_name) {

  limma::plotMA(marray_lmobject,
                coef   = input_list,
                status = geneset_results_matrix,
                values = de_call_names,
                cex    = 0.35, col = c("red", "green"),
                main   = paste0(input_list_name), legend = FALSE)
  abline(h = 0, col = "darkgrey")
  legend("bottomright", cex = 0.5, legend=(c(paste0("up in ", gene_set_list_name), paste0("down in ", gene_set_list_name))), bty =  "n", pch=c(16, 16), col = c("red",
                                                                                                                                                          "green"), text.col = c("black", "black"))
}

MAplotGeneSet3Limma <- function(marray_lmobject, geneset_results_matrix, de_call_names, input_list, input_list_name, gene_set_list_name){ 
  limma::plotMA(marray_lmobject, coef=input_list, cex=0.75, status=geneset_results_matrix, values=de_call_names, col=c("#FF0000", "#00FF00", "#507DC3"), main=paste0(input_list_name), legend=F)
  abline(h=0,col="darkgrey")
  legend("bottomright", cex = 0.5, legend=(c(paste0("up in ", gene_set_list_name), c(paste0("unchanged in ", gene_set_list_name)), paste0("down in ", gene_set_list_name))), bty =  "n", pch=c(16, 16), col = c("#FF0000", "#00FF00", "#507DC3"), text.col = c("black", "black", "black"))
}
