# Simulated customer satisfaction survey used as example dataset.
# 30 (fictitious) hotels, 40-250 guests each, 7 items + an overall score,
# all rated on a 0-10 scale. Run from the package root.
set.seed(20261021)
n_hotels <- 30
items <- c("Staff", "Cleanliness", "Comfort", "Location", "Price", "Food", "WiFi")
hotel_q <- rnorm(n_hotels, 0, 1)                       # overall quality
hotel_item <- matrix(rnorm(n_hotels * length(items), 0, 0.8), n_hotels)  # specific strengths
item_mean <- c(7.5, 7.2, 6.8, 7.8, 5.9, 6.4, 6.0)
item_load <- c(1.1, 1.0, 0.9, 0.4, 0.7, 0.9, 0.6)
polar <- runif(n_hotels, 0.6, 1.8)                     # heterogeneity of guests
rows <- list()
for (h in seq_len(n_hotels)) {
  ng <- sample(40:250, 1)
  theta <- hotel_q[h] + rnorm(ng, 0, polar[h])
  sc <- sapply(seq_along(items), function(j) {
    z <- item_mean[j] + item_load[j] * theta + hotel_item[h, j] + rnorm(ng, 0, 1.2)
    pmin(pmax(round(z), 0), 10)
  })
  colnames(sc) <- items
  overall <- pmin(pmax(round(0.35 * sc[, "Staff"] + 0.25 * sc[, "Cleanliness"] +
                               0.2 * sc[, "Comfort"] + 0.2 * sc[, "Price"] + rnorm(ng, 0, 0.8)), 0), 10)
  rows[[h]] <- data.frame(hotel = sprintf("Hotel_%02d", h), sc, Overall = overall)
}
raw <- do.call(rbind, rows)
devtools::load_all(".")
Satisfaction <- data2MatD(raw, unit = "hotel", vars = c(items, "Overall"))
Satisfaction_raw <- raw
save(Satisfaction, file = "data/Satisfaction.rda", compress = "xz")
save(Satisfaction_raw, file = "data/Satisfaction_raw.rda", compress = "xz")
