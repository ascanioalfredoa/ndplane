clamp.env2 <- function (model, env) 
{
    if (!inherits(model, "enmtools.model") & !inherits(model, 
                                                       "data.frame")) {
        stop("Argument 'model' must contain an enmtools.model object or analysis.df data frame!")
    }
    if (!inherits(env, c("SpatRaster"))) {
        stop("Argument env must be SpatRaster object!")
    }
    if (inherits(model, "enmtools.model")) {
        df <- model$analysis.df
    }
    else {
        df <- model
    }
    keep.env <- names(env)[names(env) %in% colnames(df)]
    clamped.env <- env[[keep.env]]
    for (i in 1:length(names(clamped.env))) {
        thismin <- min(df[, i], na.rm = T)
        thismax <- max(df[, i], na.rm = T)
        clamped.env[[i]] <- clamp(clamped.env[[i]], lower = thismin, upper = thismax, values = TRUE)
    }
    return(clamped.env)
}

environment(clamp.env2) <- asNamespace('ENMTools')
assignInNamespace("clamp.env", clamp.env2, ns = "ENMTools")

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

marginal.plots2 <- function (model, env, layer, standardize = TRUE, verbose = FALSE) 
{
    if (!layer %in% names(env)) {
        stop(paste("Couldn't find layer named", layer, "in environmental rasters!"))
    }
    if (inherits(model, c("enmtools.bc", "enmtools.dm"))) {
        points <- model$analysis.df[, 1:2]
    }
    else {
        points <- model$analysis.df[model$analysis.df$presence == 
                                        1, 1:2]
    }
    if (length(names(env)) == 1) {
        oldname <- names(env)
        env <- c(env, env)
        env[[2]][!is.na(env[[2]])] <- 0
        names(env) <- c(oldname, "dummyvar")
    }
    
    minmax <- terra::minmax(env)
    if (any(is.na(minmax))) {
        env <- terra::setMinMax(env)
        message("\n\nSetting min and max for environment layers...\n\n")
    }
    names <- layer
    minmax <- terra::minmax(env[[layer]])
    plot.df <- seq(minmax[1, ], minmax[2, ], length = 1000)
    for (i in names(env)) {
        if (i != layer) {
            layer.values <- terra::extract(env[[i]], points, 
                                           ID = FALSE)
            plot.df <- cbind(plot.df, rep(mean(unlist(layer.values), 
                                               na.rm = TRUE), 1000))
            names <- c(names, i)
        }
    }
    if (standardize == TRUE) {
        plot.df <- data.frame(plot.df)
    }
    colnames(plot.df) <- names
    minmax <- terra::minmax(env[[layer]])
    if (inherits(model, c("enmtools.bc", "enmtools.dm"))) {
        pres.env <- unlist(terra::extract(env[[layer]], model$analysis.df[, 
                                                                          1:2], ID = FALSE))
        pres.dens <- density(pres.env, from = minmax[1, ], to = minmax[2, 
        ], n = 1000, na.rm = TRUE)$y
        pres.dens <- pres.dens/max(pres.dens)
        bg.env <- unlist(terra::spatSample(env[[layer]], size = 10000, 
                                           na.rm = TRUE))
        bg.dens <- density(bg.env, from = minmax[1, ], to = minmax[2, 
        ], n = 1000, na.rm = TRUE)$y
        bg.dens <- bg.dens/max(bg.dens)
    }
    else {
        pres.env <- unlist(terra::extract(env[[layer]], model$analysis.df[model$analysis.df$presence == 
                                                                              1, c(1, 2)], ID = FALSE))
        pres.dens <- density(pres.env, from = minmax[1, ], to = minmax[2, 
        ], n = 1000, na.rm = TRUE)$y
        pres.dens <- pres.dens/max(pres.dens)
        bg.env <- unlist(terra::extract(env[[layer]], model$analysis.df[model$analysis.df$presence == 
                                                                            0, c(1, 2)], ID = FALSE))
        bg.dens <- density(bg.env, from = minmax[1, ], to = minmax[2, 
        ], n = 1000, na.rm = TRUE)$y
        bg.dens <- bg.dens/max(bg.dens)
    }
    if (inherits(model$model, what = "DistModel")) {
        if (verbose) {
            pred <- predict(model$model, x = plot.df, type = "response")
        }
        else {
            invisible(capture.output(pred <- predict(model$model, 
                                                     x = plot.df, type = "response")))
        }
    }
    else {
        if (inherits(model$model, "ranger")) {
            pred <- predict(model$model, data = plot.df, type = "response")$predictions[, 
                                                                                        2, drop = TRUE]
        }
        else {
            pred <- predict(model$model, newdata = plot.df, 
                            type = "response")
        }
    }
    pred <- pred/max(pred)
    plot.df.long <- data.frame(layer = c(plot.df[, layer], plot.df[, 
                                                                   layer], plot.df[, layer]), value = c(pred, pres.dens, 
                                                                                                        bg.dens), source = c(rep("Suitability", 1000), rep("Presence", 
                                                                                                                                                           1000), rep("Background", 1000)))
    response.plot <- ggplot(data = plot.df.long, aes(x = layer, 
                                                     y = value)) + geom_line(aes(colour = source, linetype = source)) + 
        xlab(layer) + ylab("Value") + theme_bw() + scale_color_manual(values = c("green4", 
                                                                                 "red", "blue")) + scale_linetype_manual(values = c("dashed", 
                                                                                                                                    "twodash", "solid")) + theme(plot.title = element_text(hjust = 0.5)) + 
        theme(plot.title = element_text(hjust = 0.5)) + theme(legend.title = element_blank())
    return(response.plot)
}

environment(marginal.plots2) <- asNamespace('ENMTools')
assignInNamespace("marginal.plots", marginal.plots2, ns = "ENMTools")
