Boxplot <- function(result, metric = "Balanced Accuracy", x = "`Classifier Name`", fill = "dataset"){
    
    #! Add checks
    
    if(!is.list(result)) result = list(result)
    
    ch <- lapply(result, function(w){
        v <- w@characteristics$value
        names(v) <- w@characteristics$characteristic
        v
    }
    )
    
    ch <- do.call("rbind", ch) |>
        as.data.frame()
    
    ch$dataset[ch$multiViewMethod != "none"] <- paste(ch$multiViewMethod, ch$dataset, sep = " - ")[ch$multiViewMethod != "none"]
    
    matrix_df = lapply(result, ClassifyR::calcCVperformance,metric) |> 
        lapply(ClassifyR::performance) |>
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
