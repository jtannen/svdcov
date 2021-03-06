setClass(
  "SVDParams",
  slots=c(
    row_cov="matrix",
    row_scores="data.frame",
    col_scores="data.frame",
    svd_d = "numeric",
    log="logical"
  )
)

SVDParams <- function(
  row_cov,
  row_scores,
  col_scores,
  svd_d,
  log=FALSE
){

  if(is.null(row.names(row_cov))) {
    row.names(row_cov) <- row_scores$row
  } else if(any(row.names(row_cov) != row_scores$row)) {
    bad_row <- which(row.names(row_cov) != row_scores$row)[1]
    stop(
      sprintf(
        "Invalid row.names for SVDParams. e.g. row %s: %s != %s",
        row.names(row_cov)[bad_row],
        row_scores$row[bad_row]
      )
    )
  }
  new(
    "SVDParams",
    row_cov=row_cov,
    row_scores=row_scores,
    col_scores=col_scores,
    svd_d=svd_d,
    log=log
  )
}
