
# Base for all graphs
common_theme <- function() { 
  theme_base() +
    theme(plot.title = element_text(face = "bold",
                                         size = rel(1.2), hjust = 0.5),
     panel.background = element_rect(colour = NA),
     plot.background = element_rect(colour = NA),
     panel.border = element_rect(colour = NA),
     axis.title = element_text(face = "bold",size = rel(1)),
     axis.title.y = element_text(angle=90,vjust =2),
     axis.title.x = element_text(vjust = -0.2),
     axis.text = element_text(size=rel(1.1)), 
     axis.line = element_line(colour="black"),
     axis.ticks = element_line(),
     panel.grid.major = element_blank(),
     panel.grid.minor = element_blank(),
     legend.key = element_rect(colour = NA),
     legend.position = "right",
     legend.key.size= unit(0.6, "cm"),
     legend.margin = unit(0, "cm"),
     legend.title = element_text(face="italic", size=rel(1)),
     plot.margin=unit(c(10,5,5,5),"mm"),
     strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
     strip.text = element_text(face="bold") 
    )
}



########################### Standard Color Palette ################################
# Followed this guy: https://drsimonj.svbtle.com/creating-corporate-colour-palettes-for-ggplot2
jens_colors <- c(
    blue = "#0073c2", 
    yellow = "#efc000", 
    grey = "#868686",
	salmon = "#cd534c", 
    cyan = "#7aa6dc", 
    dark_blue = "#003c67", 
    dark_yellow = "#8f7700", 
    dark_grey = "#3b3b3b",
	dark_salmon = "#a73030", 
    grey_blue = "#4a6990")
)

jens_colors_extract <- function(...){
    cols <- c(...)
    if (is.null(cols)) return(jens_colors)

    jens_colors[cols]
}

jens_palettes <- list(
    main = jens_colors_extract()
)

jens_pal <- function(palette="main", reverse=FALSE, ...){
    pal <- jens_palettes[[palette]]
    if (reverse) pal <- rev(pal)
    colorRampPalette(pal, ...)
}

scale_color_jens <- function(palette="main", discrete=TRUE,
  reverse=FALSE, ...) {
      pal <- jens_pal(palette=palette, reverse=reverse)

      if (discrete){
          discrete_scale("color", paste0("jens_", palette), palette=pal, ...)
      }
      else {
          scale_color_gradient(colors=pal(256), ...)
      }
    }

scale_fill_jens <- function(palette="main", discrete=TRUE, 
  reverse=FALSE, ...) {
    pal <- jens_pal(palette=palette, reverse=reverse)

    if (discrete){
        discrete_scale("fill", paste0("jens_", palette), palette=pal, ...)
    }
    else {
        scale_fill_gradientn(colors=pal(256), ...)
    }

}