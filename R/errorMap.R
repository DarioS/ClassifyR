setGeneric("errorMap", function(results, ...)
{standardGeneric("errorMap")})

setMethod("errorMap", "list", 
          function(results,
                   errorColours = list(c("#FFFFFF", "#BFBFFF", "#7F7FFF", "#3F3FFF", "#0000FF"),
                                       c("#FFFFFF", "#FFBFBF", "#FF7F7F", "#FF3F3F", "#FF0000")),
                   classColours = c("blue", "red"), fontSizes = c(24, 16, 12, 12, 12),
                   mapHeight = 4, title = "Error Comparison", showLegends = TRUE, showXlabels = TRUE,
                   showXtitle = TRUE, showYlabels = TRUE, showYtitle = TRUE,
                   legendSize = grid::unit(1, "lines"), plot = TRUE)
          {
            if(!requireNamespace("ggplot2", quietly = TRUE))
              stop("The package 'ggplot2' could not be found. Please install it.")  
            if(!requireNamespace("gridExtra", quietly = TRUE))
              stop("The package 'gridExtra' could not be found. Please install it.")       
            if(!requireNamespace("gtable", quietly = TRUE))
              stop("The package 'gtable' could not be found. Please install it.")   
            
            nColours <- if(is.list(errorColours)) length(errorColours[[1]]) else length(errorColours)
            errorBinEnds <- seq(0, 1, 1/nColours)
            knownClasses <- actualClasses(results[[1]])
            
            errors <- lapply(results, function(result)
            {
              resultTable <- do.call(rbind, predictions(result))
              sampleErrors <- by(resultTable, resultTable[, "sample"],
                                 function(samplePredictions)
                                   sum(samplePredictions[, "predicted"] != result@actualClasses[samplePredictions[1, "sample"]]))
              cut(sampleErrors / table(resultTable[, "sample"]), errorBinEnds, include.lowest = TRUE)
            })
            classedErrors <- lapply(errors, function(errorSet)
            {
              errorSet <- factor(paste(knownClasses, errorSet, sep = ','),
                                 levels = c(t(outer(levels(knownClasses), levels(errorSet), paste, sep = ','))))
            })
            
            ordering <- order(knownClasses)
            knownClasses <- knownClasses[ordering]
            errors <- lapply(errors, function(resultErrors) resultErrors[ordering])
            
            plotData <- data.frame(name = factor(rep(sampleNames(results[[1]])[ordering], length(results)), levels = sampleNames(results[[1]])[ordering]),
                                   type = factor(rep(names(results), sapply(errors, length)), levels = rev(names(results))),
                                   class = rep(knownClasses, length(results)),
                                   Error = unlist(errors))
            
            originalLegends <- showLegends                       
            originalErrorColours <- errorColours
            showLegends <- FALSE
            if(is.list(errorColours))
            {
              errorColours <- unlist(errorColours)
              plotData[, "Error"] <- unlist(classedErrors)
            }                         
            
            classData <- data.frame(Class = knownClasses)
            classesPlot <- ggplot2::ggplot(classData, ggplot2::aes(1:length(knownClasses), factor(1)), environment = environment()) +
              ggplot2::scale_fill_manual(values = classColours) + ggplot2::geom_tile(ggplot2::aes(fill = Class, height = 10)) +
              ggplot2::scale_x_discrete(expand = c(0, 0), breaks = NULL, limits = c(1, length(knownClasses))) +
              ggplot2::scale_y_discrete(expand = c(0, 0), breaks = NULL) +
              ggplot2::labs(x = '', y = '') + ggplot2::theme(plot.margin = grid::unit(c(0.2, 0, -1, 0), "lines"),
                                                             legend.title = ggplot2::element_text(size = fontSizes[4]),
                                                             legend.text = ggplot2::element_text(size = fontSizes[5]),
                                                             legend.position = ifelse(showLegends, "right", "none"),
                                                             legend.key.size = legendSize)
            
            xLabel <- if(showXtitle == TRUE) "Sample Name" else NULL
            yLabel <- if(showYtitle == TRUE) "Result" else NULL                                                             
            errorPlot <- ggplot2::ggplot(plotData, ggplot2::aes(name, type)) + ggplot2::geom_tile(ggplot2::aes(fill = Error)) +
              ggplot2::scale_fill_manual(values = errorColours, drop = FALSE) + ggplot2::scale_x_discrete(expand = c(0, 0)) +
              ggplot2::scale_y_discrete(expand = c(0, 0)) + ggplot2::theme_bw() +
              ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                             axis.text.x = if(showXlabels == TRUE) ggplot2::element_text(angle = 45, hjust = 1, size = fontSizes[3]) else ggplot2::element_blank(),
                             axis.text.y = if(showYlabels == TRUE) ggplot2::element_text(size = fontSizes[3]) else ggplot2::element_blank(),
                             axis.title.x = ggplot2::element_text(size = fontSizes[2]),
                             axis.title.y = ggplot2::element_text(size = fontSizes[2]),
                             plot.margin = grid::unit(c(0, 1, 1, 1), "lines"),
                             legend.title = ggplot2::element_text(size = fontSizes[4]),
                             legend.text = ggplot2::element_text(size = fontSizes[5]),
                             legend.position = ifelse(showLegends, "right", "none"),
                             legend.key.size = legendSize) + ggplot2::labs(x = xLabel, y = yLabel)
            
            classGrob <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(classesPlot))
            errorGrob <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(errorPlot))
            commonWidth <- grid::unit.pmax(classGrob[["widths"]], errorGrob[["widths"]])
            classGrob[["widths"]] <- commonWidth
            errorGrob[["widths"]] <- commonWidth
            
            if(originalLegends == TRUE)
            {
              showLegends <- TRUE
              errorColours <- originalErrorColours
              
              classesPlot <- ggplot2::ggplot(classData, ggplot2::aes(1:length(knownClasses), factor(1)), environment = environment()) +
                ggplot2::scale_fill_manual(values = classColours) + ggplot2::geom_tile(ggplot2::aes(fill = Class, height = 10)) +
                ggplot2::scale_x_discrete(expand = c(0, 0), breaks = NULL, limits = c(1, length(knownClasses))) +
                ggplot2::scale_y_discrete(expand = c(0, 0), breaks = NULL) +
                ggplot2::labs(x = '', y = '') + ggplot2::theme(plot.margin = grid::unit(c(0.2, 0, -1, 0), "lines"),
                                                               legend.title = ggplot2::element_text(size = fontSizes[4]),
                                                               legend.text = ggplot2::element_text(size = fontSizes[5]),
                                                               legend.position = ifelse(showLegends, "right", "none"),
                                                               legend.key.size = legendSize)
              
              classGrobUnused <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(classesPlot))            
              plotData[, "Error"] <- unlist(errors)
              errorPlot <- ggplot2::ggplot(plotData, ggplot2::aes(name, type)) + ggplot2::geom_tile(ggplot2::aes(fill = Error)) +
                ggplot2::scale_fill_manual(name = paste(levels(knownClasses)[1], "Error"),
                                           values = if(is.list(errorColours)) errorColours[[1]] else errorColours, drop = FALSE) + ggplot2::scale_x_discrete(expand = c(0, 0)) +
                ggplot2::scale_y_discrete(expand = c(0, 0)) + ggplot2::theme_bw() +
                ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                               axis.text.x = if(showXlabels == TRUE) ggplot2::element_text(angle = 45, hjust = 1, size = fontSizes[3]) else ggplot2::element_blank(),
                               axis.text.y = if(showYlabels == TRUE) ggplot2::element_text(size = fontSizes[3]) else ggplot2::element_blank(),
                               axis.title.x = ggplot2::element_text(size = fontSizes[2]),
                               axis.title.y = ggplot2::element_text(size = fontSizes[2]),
                               plot.margin = grid::unit(c(0, 1, 1, 1), "lines"),
                               legend.title = ggplot2::element_text(size = fontSizes[4]),
                               legend.text = ggplot2::element_text(size = fontSizes[5]),
                               legend.position = ifelse(showLegends, "right", "none"),
                               legend.key.size = legendSize) + ggplot2::labs(x = xLabel, y = yLabel)
              
              errorGrobUnused <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(errorPlot))
              commonWidth <- grid::unit.pmax(classGrobUnused[["widths"]], errorGrobUnused[["widths"]])                   
              errorGrobUnused[["widths"]] <- commonWidth  
              classGrobUnused[["widths"]] <- commonWidth  
              classLegend <- classGrobUnused[["grobs"]][[which(sapply(classGrobUnused[["grobs"]], function(grob) grob[["name"]]) == "guide-box")]]
              if(showLegends == TRUE)    
                firstLegend <- errorGrobUnused[["grobs"]][[which(sapply(errorGrobUnused[["grobs"]], function(grob) grob[["name"]]) == "guide-box")]]
              
              errorPlot <- ggplot2::ggplot(plotData, ggplot2::aes(name, type)) + ggplot2::geom_tile(ggplot2::aes(fill = Error)) +
                ggplot2::scale_fill_manual(name = paste(levels(knownClasses)[2], "Error"),
                                           values = if(is.list(errorColours)) errorColours[[2]] else errorColours, drop = FALSE) + ggplot2::scale_x_discrete(expand = c(0, 0)) +
                ggplot2::scale_y_discrete(expand = c(0, 0)) + ggplot2::theme_bw() +
                ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                               axis.text.x = if(showXlabels == TRUE) ggplot2::element_text(angle = 45, hjust = 1, size = fontSizes[3]) else ggplot2::element_blank(),
                               axis.text.y = if(showYlabels == TRUE) ggplot2::element_text(size = fontSizes[3]) else ggplot2::element_blank(),
                               axis.title.x = ggplot2::element_text(size = fontSizes[2]),
                               axis.title.y = ggplot2::element_text(size = fontSizes[2]),
                               plot.margin = grid::unit(c(0, 1, 1, 1), "lines"),
                               legend.title = ggplot2::element_text(size = fontSizes[4]),
                               legend.text = ggplot2::element_text(size = fontSizes[5]),
                               legend.position = ifelse(showLegends, "right", "none"),
                               legend.key.size = legendSize) + ggplot2::labs(x = xLabel, y = yLabel)
              
              errorGrobUnused <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(errorPlot)) 
              errorGrobUnused[["widths"]] <- commonWidth  
              secondLegend <- errorGrobUnused[["grobs"]][[which(sapply(errorGrobUnused[["grobs"]], function(grob) grob[["name"]]) == "guide-box")]]
            }
            
            if(showLegends == TRUE)
            {
              legendWidth <- sum(classLegend[["widths"]])
              legendHeight <- sum(firstLegend[["heights"]])
              widths <- unit.c(unit(1, "npc") - legendWidth, legendWidth)
              heights <- unit.c(unit(1 / (mapHeight + 1), "npc"), legendHeight)
              if(is.list(errorColours))
              {
                heights <- unit.c(heights, legendHeight)
                heights <- unit.c(heights, unit(1, "npc") - 2 * legendHeight - unit(1 / (mapHeight + 1), "npc"))
              } else
              {
                heights <- unit.c(heights, unit(1, "npc") - legendHeight - unit(1 / (mapHeight + 1), "npc"))
              }
            }
            else
            {
              widths <- unit(1, "npc")
              heights <- unit.c(unit(1 / (mapHeight + 1), "npc"), unit(1, "npc") - unit(1 / (mapHeight + 1), "npc"))
            }
            
            grobTable <- gtable::gtable(widths, heights)
            grobTable <- gtable::gtable_add_grob(grobTable, classGrob, 1, 1)
            if(showLegends == TRUE)
            {
              if(is.list(errorColours))
                grobTable <- gtable::gtable_add_grob(grobTable, errorGrob, 2, 1, 4, 1)
              else
                grobTable <- gtable::gtable_add_grob(grobTable, errorGrob, 2, 1, 3, 1)
            } else
            {
              grobTable <- gtable::gtable_add_grob(grobTable, errorGrob, 2, 1)
            }
            if(showLegends == TRUE)
            {  
              grobTable <- gtable::gtable_add_grob(grobTable, classLegend, 1, 2)
              grobTable <- gtable::gtable_add_grob(grobTable, firstLegend, 2, 2)
            }
            if(showLegends == TRUE && is.list(errorColours))
            {
              grobTable <- gtable::gtable_add_grob(grobTable, secondLegend, 3, 2)
              grobTable <- gtable::gtable_add_grob(grobTable, grob(), 4, 2)
            }
            wholePlot <- gridExtra::arrangeGrob(grobTable, main = grid::textGrob(title, gp = grid::gpar(fontsize = fontSizes[1])))
            
            if(plot == TRUE)               
              grid::grid.draw(wholePlot)
            wholePlot
          })
