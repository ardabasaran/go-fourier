package go_fourier

import (
	"errors"
)

var NumWorkers = 8

func dctWorker(rows <-chan int, jobReturns chan<- bool, signals [][]float64, forward bool) {
	for i := range rows {
		if forward {
			signals[i], _ = DCT1D(signals[i])
		} else {
			signals[i], _ = DCTInverse1D(signals[i])
		}
		jobReturns <- true
	}
}

func dct2D(signals [][]float64, forward bool) ([][]float64, error) {
	var err error
	height := len(signals)
	// check that input has at least one row
	if height == 0 {
		return make([][]float64, 0), errors.New("dct2D: Input 2d-array must have at least one row")
	}

	width := len(signals[0])
	// check that input has at least one column
	if width == 0 {
		return make([][]float64, 0), errors.New("dct2D: Input 2d-array must have at least one column")
	}

	// Create the result array
	result := make([][]float64, height)
	for i := 0; i < height; i++ {
		result[i] = make([]float64, width)
		for j := 0; j < width; j++ {
			result[i][j] = signals[i][j]
		}
	}

	// Apply DCT on rows as 1d arrays
	rows := make(chan int, height)
	jobReturns := make(chan bool, height)
	for w := 0; w < NumWorkers; w++ {
		go dctWorker(rows, jobReturns, result, forward)
	}

	// Send rows channel each row
	for i := 0; i < height; i++ {
		rows <- i
	}
	close(rows)
	// Wait on workers to complete
	for i := 0; i < height; i++ {
		<-jobReturns
	}
	close(jobReturns)

	// Transpose the array
	transpose := transposeReal(result)

	// Apply DFT on columns as 1d arrays
	columns := make(chan int, width)
	jobReturns = make(chan bool, width)

	for w := 0; w < NumWorkers; w++ {
		go dctWorker(columns, jobReturns, transpose, forward)
	}
	// Send columns channel each column
	for i := 0; i < width; i++ {
		columns <- i
	}
	close(columns)
	// Wait on workers to complete
	for i := 0; i < width; i++ {
		<-jobReturns
	}
	close(jobReturns)

	result = transposeReal(transpose)
	if err != nil {
		return result, err
	}

	return result, nil
}

func transposeReal(signals [][]float64) [][]float64 {
	width := len(signals)
	height := len(signals[0])
	result := make([][]float64, height)
	for i := 0; i < height; i++ {
		result[i] = make([]float64, width)
	}
	for i := 0; i < height; i++ {
		for j := 0; j < width; j++ {
			result[i][j] = signals[j][i]
		}
	}
	return result
}

func dftWorker(rows <-chan int, jobReturns chan<- bool, signals [][]complex128, forward bool, algorithm string) {
	for i := range rows {
		if forward {
			switch algorithm {
			case "radix2":
				DFT2Radix1D(signals[i])
			default:
				DFTNaive1D(signals[i])
			}

		} else {
			switch algorithm {
			case "radix2":
				DFTInverse2Radix1D(signals[i])
			default:
				DFTInverseNaive1D(signals[i])
			}
		}
		jobReturns <- true
	}
}

func dft2D(signals [][]complex128, forward bool, algorithm string) error {
	var err error
	height := len(signals)
	// check that input has at least one row
	if height == 0 {
		return errors.New("dft2D: Input 2d-array must have at least one row")
	}

	width := len(signals[0])
	// check that input has at least one column
	if width == 0 {
		return errors.New("dft2D: Input 2d-array must have at least one column")
	}

	// Apply DFT on rows as 1d arrays
	rows := make(chan int, height)
	jobReturns := make(chan bool, height)
	for w := 0; w < NumWorkers; w++ {
		go dftWorker(rows, jobReturns, signals, forward, algorithm)
	}

	// Send rows channel each row
	for i := 0; i < height; i++ {
		rows <- i
	}
	close(rows)
	// Wait on workers to complete
	for i := 0; i < height; i++ {
		<-jobReturns
	}
	close(jobReturns)

	// Transpose the array
	transpose := transposeComplex(signals)

	// Apply DFT on columns as 1d arrays
	columns := make(chan int, width)
	jobReturns = make(chan bool, width)

	for w := 0; w < NumWorkers; w++ {
		go dftWorker(columns, jobReturns, transpose, forward, algorithm)
	}
	// Send columns channel each column
	for i := 0; i < width; i++ {
		columns <- i
	}
	close(columns)
	// Wait on workers to complete
	for i := 0; i < width; i++ {
		<-jobReturns
	}
	close(jobReturns)

	transpose = transposeComplex(transpose)
	if err != nil {
		return err
	}
	for i := 0; i < height; i++ {
		for j := 0; j < width; j++ {
			signals[i][j] = transpose[i][j]
		}
	}
	return nil
}

func transposeComplex(signals [][]complex128) [][]complex128 {
	width := len(signals)

	height := len(signals[0])

	// create the result array
	result := make([][]complex128, height)

	for i := 0; i < height; i++ {
		result[i] = make([]complex128, width)
	}

	for i := 0; i < height; i++ {
		for j := 0; j < width; j++ {
			result[i][j] = signals[j][i]
		}
	}

	return result
}