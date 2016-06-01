\name{SlicedData-class}
\Rdversion{1.1}
\docType{class}
\alias{SlicedData-class}
\alias{SlicedData}
\alias{SlicedData-class}
\alias{[[,SlicedData-method}
\alias{[[<-,SlicedData-method}
\alias{colnames,SlicedData-method}
\alias{colnames<-,SlicedData-method}
\alias{dim,SlicedData-method}
\alias{length,SlicedData-method}
\alias{ncol,SlicedData-method}
\alias{NCOL,SlicedData-method}
\alias{nrow,SlicedData-method}
\alias{NROW,SlicedData-method}
\alias{rownames,SlicedData-method}
\alias{rownames<-,SlicedData-method}
\alias{show,SlicedData-method}
\alias{as.matrix,SlicedData-method}
\alias{rowMeans,SlicedData-method}
\alias{rowSums,SlicedData-method}
\alias{colMeans,SlicedData-method}
\alias{colSums,SlicedData-method}
\alias{summary.SlicedData}

\title{Class \code{SlicedData} for storing large matrices}
\description{
	This class is created for fast and memory efficient manipulations with large datasets presented in matrix form.
	It is used to load, store, and manipulate large datasets, e.g. genotype and gene expression matrices.
	When a dataset is loaded, it is sliced in blocks of 1,000 rows (default size).
	This allows imputing, standardizing, and performing other operations with the data with minimal memory overhead.
}
\value{
check the resutls
}
\section{Extends}{
	\code{SlicedData} is a reference classes (\code{\linkS4class{envRefClass}}).
	Its methods can change the values of the fields of the class.
}
\references{
	The package website: \url{http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/}
}
\author{
	Andrey Shabalin \email{ashabalin@vcu.edu}
}
\keyword{classes}
\section{Fields}{
	\describe{
		\item{\code{dataEnv}:}{\code{environment}. Stores the slices of the data matrix. The slices should be accessed via \code{getSlice()} and \code{setSlice()} methods. }
		\item{\code{nSlices1}:}{\code{numeric}. Number of slices. For internal use. The value should be access via \code{nSlices()} method. }
		\item{\code{rowNameSlices}:}{\code{list}. Slices of row names. }
		\item{\code{columnNames}:}{\code{character}. Column names. }
		\item{\code{fileDelimiter}:}{\code{character}. Delimiter separating values in the input file. }
		\item{\code{fileSkipColumns}:}{\code{numeric}. Number of columns with row labels in the input file. }
		\item{\code{fileSkipRows}:}{\code{numeric}. Number of rows with column labels in the input file. }
		\item{\code{fileSliceSize}:}{\code{numeric}. Maximum number of rows in a slice. }
		\item{\code{fileOmitCharacters}:}{\code{character}. Missing value (NaN) representation in the input file. }
	}
}
\section{Methods}{
	\describe{
		\item{\code{initialize(mat)}:}{ Create the object from a matrix. }
		\item{\code{nSlices()}:}{ Returns the number of slices. }
		\item{\code{nCols()}:}{ Returns the number of columns in the matrix. }
		\item{\code{nRows()}:}{ Returns the number of rows in the matrix. }
		\item{\code{Clear()}:}{ Clears the object. Removes the data slices and row and column names. }
		\item{\code{Clone()}:}{ Makes a copy of the object. Changes to the copy do not affect the source object. }
		\item{\code{CreateFromMatrix(mat)}:}{ Creates \code{SlicedData} object from a \code{\link[base]{matrix}}.}
		\item{\code{LoadFile(filename, skipRows = NULL, skipColumns = NULL,} \cr 
			\code{sliceSize = NULL, omitCharacters = NULL, delimiter = NULL, rowNamesColumn = 1)}:}{ Loads data matrix from a file. \code{filename} should be a character string. The remaining parameters specify the file format and have the same meaning as \code{file*} fields. Additional \code{rowNamesColumn} parameter specifies which of the columns of row labels to use as row names.}
		\item{\code{SaveFile(filename)}:}{ Saves the data to a file. \code{filename} should be a character string.}
		\item{\code{getSlice(sl)}:}{ Retrieves \code{sl}-th slice of the matrix. }
		\item{\code{setSlice(sl, value)}:}{ Set \code{sl}-th slice of the matrix. }
		\item{\code{ColumnSubsample(subset)}:}{ Reorders/subsets the columns according to \code{subset}. \cr
Acts as \code{M = M[ ,subset]} for a matrix \code{M}. }
		\item{\code{RowReorder(ordr)}:}{ Reorders rows according to \code{ordr}. \cr
Acts as \code{M = M[ordr, ]} for a matrix \code{M}. }
		\item{\code{RowMatrixMultiply(multiplier)}:}{ Multiply each row by the \code{multiplier}. \cr
Acts as \code{M = M \%*\% multiplier} for a matrix \code{M}. }
		\item{\code{CombineInOneSlice()}:}{ Combines all slices into one. The whole matrix can then be obtained via \code{$getSlice(1)}. }
		\item{\code{IsCombined()}:}{ Returns \code{TRUE} if the number of slices is 1 or 0. }
		\item{\code{ResliceCombined(sliceSize = -1)}:}{ Cuts the data into slices of \code{sliceSize} rows. If \code{sliceSize} is not defined, the value of \code{fileSliceSize} field is used.}
		\item{\code{GetAllRowNames()}:}{ Returns all row names in one vector. }	
		\item{\code{RowStandardizeCentered()}:}{ Set the mean of each row to zero and the sum of squares to one. }
		\item{\code{SetNanRowMean()}:}{ Impute rows with row mean. Rows full of NaN values are imputed with zeros. }
		\item{\code{RowRemoveZeroEps()}:}{ Removes rows of zeros and those that are nearly zero. }
		\item{\code{FindRow(rowname)}:}{ Finds row by name. Returns a pair of slice number an row number within the slice. If no row is found, the function returns \code{NULL}. }
		\item{\code{rowMeans(x, na.rm = FALSE, dims = 1L)}:}{Returns a vector of row means. Works as \link[base]{rowMeans} but requires \code{dims} to be equal to \code{1L}.}
		\item{\code{rowSums(x, na.rm = FALSE, dims = 1L)}:}{Returns a vector of row sums. Works as \link[base]{rowSums} but requires \code{dims} to be equal to \code{1L}.}
		\item{\code{colMeans(x, na.rm = FALSE, dims = 1L)}:}{Returns a vector of column means. Works as \link[base]{colMeans} but requires \code{dims} to be equal to \code{1L}.}
		\item{\code{colSums(x, na.rm = FALSE, dims = 1L)}:}{Returns a vector of column sums. Works as \link[base]{colSums} but requires \code{dims} to be equal to \code{1L}.}
	}
}

\usage{
# x[[i]] indexing allows easy access to individual slices.
# It is equivalent to x$GetSlice(i) and x$SetSlice(i,value)
\S4method{[[}{SlicedData}(x, i)
\S4method{[[}{SlicedData}(x, i) <- value

# The following commands work as if x was a simple matrix object
\S4method{nrow}{SlicedData}(x)
\S4method{ncol}{SlicedData}(x)
\S4method{dim}{SlicedData}(x)
\S4method{rownames}{SlicedData}(x)
\S4method{colnames}{SlicedData}(x)
\S4method{rownames}{SlicedData}(x) <- value
\S4method{colnames}{SlicedData}(x) <- value

# SlicedData object can be easily transformed into a matrix
# preserving row and column names
\S4method{as.matrix}{SlicedData}(x)

# length(x) can be used in place of x$nSlices()
# to get the number of slices in the object
\S4method{length}{SlicedData}(x)
}
\arguments{
	\item{x}{
		\code{\linkS4class{SlicedData}} object.
	}
	\item{i}{
		Number of a slice.
	}
	\item{value}{
		New content for the slice / new row or column names.
	}
}


\examples{

# Create a SlicedData variable
}