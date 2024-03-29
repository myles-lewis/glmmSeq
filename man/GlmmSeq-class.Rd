% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClass.R
\docType{class}
\name{GlmmSeq-class}
\alias{GlmmSeq-class}
\title{An S4 class to define the glmmSeq output}
\description{
An S4 class to define the glmmSeq output
}
\section{Slots}{

\describe{
\item{\code{info}}{List including the matched call, dispersions, offset, designMatrix}

\item{\code{formula}}{The model formula}

\item{\code{stats}}{Statistics from fitted models}

\item{\code{predict}}{Predicted values}

\item{\code{reduced}}{Optional reduced formula for LRT}

\item{\code{countdata}}{The input expression data with count data in rows}

\item{\code{metadata}}{The input metadata}

\item{\code{modelData}}{Model data for predictions}

\item{\code{optInfo}}{Information on whether the model was singular or converged}

\item{\code{errors}}{Any errors}

\item{\code{vars}}{List of variables stored from the original call, including the
\code{id} variable (by default automatically identified from the random effect
term in the model) and \code{removeSingles} argument}
}}

