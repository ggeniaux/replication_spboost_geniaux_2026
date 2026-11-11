suppressMessages(library(sf))
prep <- readRDS("lucas_prep.rds")
crs_h <- st_crs(prep$meta$crs)
cty <- st_read("tiger_raw/cb_county", quiet = TRUE)
lucas <- cty[cty$STATEFP == "39" & cty$COUNTYFP == "095", ]
wat <- st_read("tiger_raw/water", quiet = TRUE)
wat <- wat[wat$AWATER > 2e5, ]              # large water bodies only
to_km <- function(x) {
  g <- st_transform(st_geometry(x), crs_h)
  st_set_crs(g / 1000, NA)
}
bm <- list(boundary = to_km(lucas), water = to_km(wat))
saveRDS(bm, "lucas_basemap.rds")
cat("boundary bbox (km):", round(st_bbox(bm$boundary), 1), "\n")
cat("water polygons:", length(bm$water), "\n")
