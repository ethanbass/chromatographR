# These functions are copied with slight modification from the R plotly package.
# Copyright (c) 2017, Plotly, Inc.
#   
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
#                                                             "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#   
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# Sys.setenv(NOT_CRAN = "true")

#' @noRd
visual_testing <- grepl("true", Sys.getenv("VISUAL_TESTS"), fixed = TRUE)

# message("Visual testing is ", if (!visual_testing) "not ", "enabled.")
# start up the orca image server
imageServer <- if (visual_testing) {
  # https://github.com/plotly/plotly.R/issues/2179
  if (!reticulate::py_module_available("kaleido")){
    reticulate::py_require("kaleido")
  }
  if (!reticulate::py_module_available("plotly")){
    reticulate::py_require("plotly")
  }
  reticulate::py_run_string("import sys") 
  plotly::kaleido() 
} else {
  list(transform = function(...) stop("Visual testing is disabled!"))
}

#' @noRd
expect_doppelganger_plotly <- function(name, p, ...) {
  testthat::local_edition(3)
  
  name <- str_standardise(name)
  file <- paste0(name, ".svg")
  path <- tempfile(file, fileext = ".svg")
  testthat::announce_snapshot_file(path, name = file)
  
  if (!visual_testing) {
    return(invisible(NULL))
  }
  
  # some plots have random characteristics, so make sure we always have the same seed,
  # otherwise comparing svg produces false positives
  set.seed(555)
  
  write_plotly_svg(p, path)
  testthat::expect_snapshot_file(
    path = path, name = file, cran = FALSE,
    compare = function(old, new) {
      compare_file_text(old, new) || identical(rsvg::rsvg_png(old), rsvg::rsvg_png(new))
    }
  )
}

#' Run visual test and return 'built' data/layout
#' @noRd
expect_doppelganger_built <- function(p, name, ...) {
  vdiffr::expect_doppelganger(p, name, ...)
  plotly::plotly_build(p)$x[c("data", "layout")]
}


#' Define logic for writing svg
#' @noRd
write_plotly_svg <- function(p, file) {
  # before exporting, specify trace[i].uid so resulting svg is deterministic
  # https://github.com/plotly/orca/issues/133
  p <- plotly::plotly_build(p)
  uid_data <- paste0("-vdiffr-plotly-", seq_along(p$x$data))
  p$x$data <- Map(function(tr, id) { tr$uid <- id; tr }, p$x$data, uid_data)
  
  # write svg to disk
  owd <- setwd(dirname(file))
  on.exit(setwd(owd))
  imageServer$transform(p, file = basename(file), width = 640, height = 480)
  
  # strip out non-deterministic fullLayout.uid
  # TODO: if and when plotly provides an API to pre-specify, use it!
  svg_txt <- readLines(file, warn = FALSE)
  strextract <- function(str, pattern) regmatches(str, regexpr(pattern, str))
  def <- strextract(svg_txt, 'defs id=\\"defs-[[:alnum:]]+\\"')
  uid <- sub("defs-", "", strextract(def, "defs-[[:alnum:]]+"))
  svg_txt <- gsub(uid, "", svg_txt, fixed = TRUE)
  writeLines(svg_txt, file)
}

#' copied from vdiffr
#' @noRd
str_standardise <- function(s, sep = "-") {
  stopifnot(is.character(s) && length(s) == 1)
  s <- gsub("[^a-z0-9]", sep, tolower(s))
  s <- gsub(paste0(sep, sep, "+"), sep, s)
  s <- gsub(paste0("^", sep, "|", sep, "$"), "", s)
  s
}
