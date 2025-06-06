% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_emuFit.R
\name{plot.emuFit}
\alias{plot.emuFit}
\title{Plotting function}
\usage{
\method{plot}{emuFit}(
  x,
  plot_key = NULL,
  title = NULL,
  taxon_names = NULL,
  display_taxon_names = TRUE,
  data_only = FALSE,
  ...
)
}
\arguments{
\item{x}{Output from emuFit()}

\item{plot_key}{(Optional) Default \code{NULL}. List of named vectors containing names in the "covariate" column of the \code{coef} output of the radEmu model object. If you wish for multiple covariate values to be plotted on the same plot, then those variables should be included in the same named vector. By default, each column of the design matrix receives its own plot.}

\item{title}{(Optional). Default \code{NULL}. Character string. The main title for the graphic.}

\item{taxon_names}{(Optional). Default \code{NULL}. Data frame. If \code{NULL}, keep taxon names as listed in radEmu model. Otherwise, users can input a data frame with two columns: one labelled "category" with the same levels as in the radEmu output and another labelled "cat_small" with the preferred labels.}

\item{display_taxon_names}{(Optional). Default \code{TRUE}. Boolean. If \code{FALSE}, remove sample names from the plot.}

\item{data_only}{(Optional). Default \code{FALSE}. Boolean. If \code{TRUE}, only returns data frame.}

\item{...}{There are no optional parameters at this time.}
}
\value{
Object of class \code{ggplot}. Plot of \code{radEmu} model fit with 95\% confidence intervals.
}
\description{
Plotting function
}
\examples{
data(wirbel_sample_small)
data(wirbel_otu_small)

ch_fit <- emuFit(formula = ~ Group + Study, 
                 data = wirbel_sample_small,
                 Y = wirbel_otu_small,
                 run_score_tests = FALSE,
                 tol = 0.01)

plot_key <- list(p1 = c("Control" = "GroupCTR"),
                 p2 = c("FR-Control" = "StudyFR-CRC",
                        "US-Control" = "StudyUS-CRC"))

out <- plot(x = ch_fit,
            plot_key = plot_key,
            display_taxon_names = FALSE)

out$plots$p1
out$plots$p2
}
