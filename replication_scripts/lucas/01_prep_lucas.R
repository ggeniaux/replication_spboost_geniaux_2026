## Prepares the Lucas County data (spData::house) and the a priori menu of W.
## Usage: Rscript 01_prep_lucas.R          (from replication_scripts/lucas)
## Output: lucas_prep.rds  (list: d, Wlist = {knn4, knn10, dband}, meta)

suppressMessages({
  library(spData); library(spdep); library(sf); library(Matrix)
})
data(house)

cc <- sp::coordinates(house)
d  <- house@data

## Toledo CBD (long/lat WGS84) projected into the data CRS
cbd_ll <- sf::st_sfc(sf::st_point(c(-83.5379, 41.6529)), crs = 4326)
crs_house <- sf::st_crs(sp::proj4string(house))
cbd <- sf::st_coordinates(sf::st_transform(cbd_ll, crs_house))
dcbd_km <- sqrt((cc[, 1] - cbd[1])^2 + (cc[, 2] - cbd[2])^2) / 1000

dat <- data.frame(
  lprice     = log(d$price),
  yrbuilt    = as.numeric(d$yrbuilt),
  TLA        = as.numeric(d$TLA),
  lotsize    = as.numeric(d$lotsize),
  frontage   = as.numeric(d$frontage),
  depth      = as.numeric(d$depth),
  garagesqft = as.numeric(d$garagesqft),
  rooms      = as.numeric(d$rooms),
  beds       = as.numeric(d$beds),
  baths      = as.numeric(d$baths),
  halfbaths  = as.numeric(d$halfbaths),
  ## rare levels (< 300 sales) pooled: avoids crossfit folds
  ## lacking the level in training (and unstable coefficients)
  wall       = factor(ifelse(
                 d$wall == "stone", "brick",
                 ifelse(d$wall %in% c("ccbtile", "stucdrvt"), "other",
                        as.character(d$wall)))),
  garage     = factor(ifelse(d$garage %in% c("basement", "carport"),
                             "other", as.character(d$garage))),
  stories    = factor(ifelse(d$stories %in% c("two+half", "three"),
                             "two", as.character(d$stories))),
  syear      = d$syear,
  dcbd       = dcbd_km,
  x          = cc[, 1] / 1000,
  y          = cc[, 2] / 1000)
stopifnot(!anyNA(dat$lprice), nrow(dat) == 25357L)

## ---- W menu fixed a priori ----------------------------------------------
lw2sparse <- function(lw) as(spatialreg::as_dgRMatrix_listw(lw), "CsparseMatrix")
knn_W <- function(k) {
  nb <- spdep::knn2nb(spdep::knearneigh(cc, k = k))
  lw2sparse(spdep::nb2listw(nb, style = "W"))
}
W4  <- knn_W(4L)
W10 <- knn_W(10L)

## distance band: smallest radius leaving no observation isolated
d1   <- unlist(spdep::nbdists(spdep::knn2nb(spdep::knearneigh(cc, k = 1L)), cc))
dmax <- max(d1)
## negative lower bound: includes pairs at zero distance
## (repeat sales of the same parcel, duplicated coordinates)
nb_db <- spdep::dnearneigh(cc, -1, dmax * 1.001)  # margin: floating-point bound
card  <- spdep::card(nb_db)
Wdb <- lw2sparse(spdep::nb2listw(nb_db, style = "W", zero.policy = FALSE))

meta <- list(n = nrow(dat), dmax_m = dmax,
             deg_dband = summary(card),
             nnz = c(knn4 = Matrix::nnzero(W4), knn10 = Matrix::nnzero(W10),
                     dband = Matrix::nnzero(Wdb)),
             cbd_xy_km = as.numeric(cbd) / 1000,
             crs = sp::proj4string(house))

saveRDS(list(d = dat, Wlist = list(knn4 = W4, knn10 = W10, dband = Wdb),
             meta = meta), "lucas_prep.rds")
cat(sprintf("n=%d | dmax=%.0f m | dband degree: med=%d max=%d | nnz: %s\n",
            meta$n, dmax, median(card), max(card),
            paste(names(meta$nnz), meta$nnz, sep = "=", collapse = " ")))
