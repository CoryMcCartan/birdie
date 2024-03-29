
#' A pseudo-voterfile
#'
#' A dataset containing 5,000 fake voter records.
#' Created by randomizing a subset of the North Carolina voter file.
#' Turnout records are completely randomly generated.
#'
#' @format A data frame with 5,000 rows and 4 records:
#' \describe{
#'   \item{last_name}{Voter's last name}
#'   \item{zip}{5-digit ZIP code. May be NA}
#'   \item{race}{One of "white", "black", "hisp", "asian", "aian", or "other"}
#'   \item{turnout}{1 if the voter voted in the most recent election, 0 otherwise}
#' }
#'
#' @source \url{https://www.ncsbe.gov/results-data/voter-registration-data}
#'
#' @concept misc
#' @examples
#' data(pseudo_vf)
#' print(pseudo_vf)
"pseudo_vf"
