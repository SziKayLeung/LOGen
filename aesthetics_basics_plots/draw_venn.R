
suppressMessages(library(VennDiagram))
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

twovenndiagrams <- function(set1, set2, name1, name2){
  p <- venn.diagram(x = list(set1,set2), 
                    category.names = c(name1,name2),filename = NULL, output=TRUE, lwd = 0.2,lty = 'blank', 
                    fill = c("#B3E2CD", "#FDCDAC"), main = "\n", cex = 1,fontface = "bold",fontfamily = "ArialMT",#
                    cat.cex = 1,  cat.default.pos = "outer",  cat.pos = c(-120, 120), cat.dist = c(0.055, 0.055),  cat.fontfamily = "ArialMT",  #rotation = 1,   main = "\n\n\n\n"
                    print.mode = "raw")
  return(p)
}
