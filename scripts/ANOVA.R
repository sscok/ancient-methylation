#!/usr/bin/env Rscript
# Description: This function is actually a built-in function that captures 
# the analysis of variance summary output in our script below. 
# However, in our system it needed to be defined.
# inputs: 
#	- expr
#     - output file name: string
#     - append: boolean
#     - collapse: string
# output:
#	- file output
captureOutput <- function(expr, file=NULL, append=FALSE, collapse=NULL, envir=parent.frame()) {
      # Argument 'file':
      # Default is to capture via a raw connection
      if (is.null(file)) file <- raw(0L)

      # It is still possible to capture via a string
      if (identical(file, character(0L))) file <- NULL

      # How to capture output?
      if (is.raw(file)) {
            # Via a temporary raw connection? [MUCH FASTER]
            res <- eval({
                        file <- rawConnection(raw(0L), open="w")
                        on.exit({
                              if (!is.null(file)) close(file)
                              })

                        capture.output(expr, file=file)

                        res <- rawConnectionValue(file)
                        close(file)
                        file <- NULL; # Not needed anymore

                        # Convert to character
                        res <- rawToChar(res)

                        res
                  }, envir=envir, enclos=envir)
      } else {
            # Backward compatibility, i.e. capture to file
            res <- eval({
                        capture.output(expr, file=file, append=append)
                  }, envir=envir, enclos=envir)

            return(invisible(res))
      }

      ## At this point 'res' is a single character string if captured
      ## to a raw or file connection, whereas if captured to say
      ## "text" connection, then it is a character vector with elements
      ## split by '\n' newlines.
      ## In order to emulate capture.output() behavior as far as possible,
      ## we will split by '\n'.
      res <- unlist(strsplit(res, split="\n", fixed=TRUE), use.names=FALSE)

      ## Merge back using the collapse string?
      if (!is.null(collapse)) res <- paste(res, collapse=collapse)

      res
} # captureOutput()


args = commandArgs(trailingOnly=TRUE)
input=read.table(args[1], head=F)		# path to replicated file for all individuals

genes=c()
for (c_gene in unique(input$V6)){ 	# per gene
 gene = input[input$V6 == c_gene,] 	# get part related to the gene in concern
 temp=tryCatch(summary(aov(gene$V3 ~ gene$V5 + gene$V7 + gene$V8 + Error(gene$V4))), error = function(e) "NaN")	# perform ANOVA
 captureOutput(temp, file = args[2], append=TRUE)     # capture the result to a file
 genes=append(genes,c_gene)
}
write.table(genes,args[3],col.names=F, row.names=F, quote=F) # get the corresponding genes that had ANOVA values to another file.
# process the files into a desired format merge the genes with the p-values to get a p-value table for all factors
