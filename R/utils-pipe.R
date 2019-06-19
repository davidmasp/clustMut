#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL



ascii_message <- function(version = packageVersion("clustMut")) {

nchar_version = nchar(as.character(version))

nchar_rest = 46 - 12 - nchar_version -1
rest_string = paste(rep(' ',nchar_rest),collapse = "")
vector = c("+--------------------------------------------+",
  "|                                            |",
  glue::glue("|  ClustMut {version}{rest_string}|"),
  "|                                            |",
  "|          +  +  +         +          +      |",
  "|          |  |  |         |          |      |",
  "|   -------+--+--+---------+----------+---   |",
  "|                                            |",
  "+--------------------------------------------+")

  for (i in vector){
    print(i)
  }


}
