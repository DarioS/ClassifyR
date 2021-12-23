Boxplot = function(result, metric = "Balanced Accuracy", x = "`Classifier Name`", fill = "dataset"){
    
    ch <- lapply(result, function(w){
        v <- w@characteristics$value
        names(v) <- w@characteristics$characteristic
        v
    }
    )
    
    ch <- do.call("rbind", ch) |>
        as.data.frame()
    
    matrix_df = lapply(result, ClassifyR::calcCVperformance,metric) |> 
        lapply(ClassifyR::performance) |>
        purrr::flatten() |> 
        rlist::list.rbind() 
    
    
    colnames(matrix_df) = "perfMetric"
    matrix_df <- cbind(ch, matrix_df)
    
        ggplot2::ggplot(matrix_df, ggplot2::aes_string(x = x  , y = "perfMetric", fill = fill)) + 
        ggplot2::geom_boxplot() + 
        ggplot2::theme_minimal() + 
        ggplot2::labs(y = metric)
}
