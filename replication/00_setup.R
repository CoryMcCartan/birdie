suppressMessages({
    library(here)
    devtools::load_all(here())
    library(tidyverse)
    library(scales)
    library(wacolors)
    library(patchwork)
    library(geomtextpath)
})

races = c(white="White", black="Black", hisp="Hispanic", asian="Asian",
          aian="Native", other="Other")

theme_paper = function() {
    theme_bw(base_family="Times", base_size=10) +
        theme(plot.background=element_blank())
}
