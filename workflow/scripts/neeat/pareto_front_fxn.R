pareto_front <- function(df, col1="false_pos",col2="false_neg") {

    res <- data.frame(as.list(df[FALSE,]))

    for (i in 1:nrow(df)) {
        front <- TRUE
        for (j in 1:nrow(df)) {
            if (i==j)
                next
            if (df[i,col1]>df[j,col1] && df[i,col2]>df[j,col2]) {
                front <- FALSE
                break
            }
        }

        if (front)
            res <- rbind(res,df[i,])
    }

    res
}

