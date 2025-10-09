Melting_x <- function(Res, x= 15, metadatas_selected_col)
{
  Res_melt <- reshape2::melt(Res[,c(1,x:dim(Res)[2])], id = "clusters.id", variable.name = "Run_Barcod", value.name = "Nb.reads")
  Res_melt <- Res_melt[Res_melt$Nb.reads != 0,]
  Res_melt <- merge(Res_melt,Res[,c(1:(x-1))], by = "clusters.id", all = T)
  Res_melt <- merge(Res_melt, metadatas[,metadatas_selected_col], by = "Run_Barcod", all = T)
  # Remove samples that haven't any result
  Res_melt <- Res_melt[!is.na(Res_melt$clusters.id),]
  # Creating Taxon column
  Res_melt[which(is.na(Res_melt$Family)),c("Family","Genus","Species")] <- "unknown"
  Res_melt$Taxon <- ifelse(is.na(Res_melt$Genus), Res_melt$Family,paste(Res_melt$Genus, Res_melt$Species))
  #unique(Res_melt$Taxon)
  return(Res_melt)
}

GAM_model_plot <- function(data, k, y, ylab)
{
  p <- ggplot(data, aes_string(x = "Distance.to.Reef..Nautic.Miles.", y = y, color = "Sampling.Site")) +
    geom_point(size = 2.5, alpha = 0.8) +
    geom_smooth(method = "gam", formula = y ~ s(x, k = k), se = TRUE,
                color = "black", fill = "grey70", alpha = 0.4, size = 1.2) +
    scale_color_brewer(palette = "Dark2") +
    scale_x_continuous(breaks = c(0, 1, 2, 5, 10)) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = c(0.8, 0.8),
      legend.background = element_rect(fill = "white", color = "grey70"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black")
    ) +
    labs(
      x = "Distance from the reef (nautical miles)",
      y = ylab,
      color = "Transect"
    ) 
  
  return(p)
}

# --- Function to repel overlapping pies (simple iterative approach) ---
repel_points_simple <- function(x, y, radius, repel_strength = 3, max_iter = 500) {
  # Minimum allowed distance between pies (average radius × factor)
  min_dist <- mean(radius) * repel_strength
  
  coords <- data.frame(X2 = x, Y2 = y)
  n <- nrow(coords)
  
  for (iter in 1:max_iter) {
    moved <- FALSE
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        dx <- coords$X2[i] - coords$X2[j]
        dy <- coords$Y2[i] - coords$Y2[j]
        dist <- sqrt(dx^2 + dy^2)
        if (dist < min_dist) {
          # Apply small repulsion shift
          shift <- (min_dist - dist) / 2
          angle <- atan2(dy, dx)
          coords$X2[i] <- coords$X2[i] + cos(angle) * shift
          coords$Y2[i] <- coords$Y2[i] + sin(angle) * shift
          coords$X2[j] <- coords$X2[j] - cos(angle) * shift
          coords$Y2[j] <- coords$Y2[j] - sin(angle) * shift
          moved <- TRUE
        }
      }
    }
    if (!moved) break
  }
  coords
}