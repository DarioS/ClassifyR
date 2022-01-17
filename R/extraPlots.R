Boxplot <- function(result, metric = "Balanced Accuracy", x = "`Classifier Name`", fill = "dataset"){
    
    #! Add checks
    
    if(!is.list(result)) result = list(result)
    
    if(length(result)>1 & sd(unlist(lapply(result,function(x)nrow(x@predictions))))>0)stop("At the moment, all results have the same number of folds and repeats")
    
    ch <- lapply(result, function(w){
        v <- w@characteristics$value
        names(v) <- w@characteristics$characteristic
        v
    }
    )
    
    ch <- dplyr::bind_rows(ch) |>
        as.data.frame()
    
    ch$dataset[ch$multiViewMethod != "none"] <- paste(ch$multiViewMethod, ch$dataset, sep = " - ")[ch$multiViewMethod != "none"]
    
    matrix_df = lapply(result, calcCVperformance,metric) |> 
        lapply(performance) |>
        purrr::flatten() |> 
        lapply(t) |>
        rlist::list.rbind() |>
        as.vector()
    
    

    matrix_df <- suppressWarnings(cbind(perfMetric = matrix_df, ch))
    
        ggplot2::ggplot(matrix_df, ggplot2::aes_string(x = x  , y = "perfMetric", fill = fill)) + 
        ggplot2::geom_boxplot() + 
        ggplot2::theme_minimal() + 
        ggplot2::labs(y = metric)
}
