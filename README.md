
# svdcov

<!-- badges: start -->
<!-- badges: end -->

The goal of `svdcov` is to calculate the covariance matrix of rectangular data 
using SVD or shrinkage methods.

## Installation

You can install the released version of svdcov from github with:

``` r
devtools::intall_github("github.io/jtannen/svdcov")
```

## Example


``` r
library(svdcov)

# Consider a block-diagonal example.
cov <- matrix(0, 10, 10)
cov[row(cov) <=5 & col(cov) <= 5] <- 0.5
cov[row(cov) > 5 & col(cov) > 5] <- 0.5
diag(cov) <- 1

data <- MASS::mvrnorm(
  n=100, 
  mu=1:10, 
  Sigma=cov
)

data <- t(data)

res <- get_svd(data, n_svd=2, col_means=rep(0, 100))

res@row_cov
get_fitted(res)

sample_from_posterior(
  res,
  obs=c(1.42, 7.34),
  obs_id=c(1, 6),
  col_mean=0,
  col_mean_sd=1
)


## Example with Philadelphia elections data
library(dplyr)
library(sf)
library(ggplot2)

data("phila_elections")
data("phila_divs")

phila_elections <- phila_elections %>%
  filter(
    !grepl("WRITE(\\s|-)IN", candidate, ignore.case=T),
    substr(warddiv, 1, 2) != "99", # overseas
    is.na(district) # city-wide races
  )

competitive_races <- phila_elections %>% 
  group_by(office, election_type, year, candidate) %>%
  summarise(votes=sum(votes), .groups="drop") %>%
  group_by(office, election_type, year) %>%
  summarise(is_competitive = sum(votes > 100) > 2) %>%
  filter(is_competitive)
  
# mat is a matrix with division rows, and candidate columns
mat <- phila_elections %>% 
  inner_join(competitive_races) %>%
  mutate(column_name = paste(candidate, office, party, election_type, year, sep="::")) %>%
  group_by(office, year, election_type, warddiv) %>%
  mutate(log_pvote = log((votes + 1) / sum(votes+1))) %>%
  ungroup() %>%
  pivot_wider(
    id_cols = "warddiv", 
    names_from="column_name", 
    values_from="log_pvote",
    values_fill=list(log_pvote=log(0.001))
  )
  
mat <- mat %>% tibble::column_to_rownames("warddiv")
mat <- as.matrix(mat)

row_ward <- substr(row.names(mat), 1, 2)
res <- get_svd(
  mat, 
  n_svd=3,
  row_groups=row_ward  ## cluster uncertainty by ward
)

ggplot(
  divs %>% left_join(res@row_scores, by=c("warddiv"="row"))
) + 
  geom_sf(aes(fill=score.1), color=NA) +
  scale_fill_gradient2() +
  theme_minimal()
  
res@col_scores %>% arrange(desc(score.1)) %>% head()

```

