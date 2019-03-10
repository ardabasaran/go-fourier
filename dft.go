package go_fourier

import (
	"errors"
	"math"
	"math/bits"
)

// DFT2Radix1D computes the discrete fourier transform of the given array in the complex number space.
// The calculation is done in place using Cooley-Tukey radix-2 algorithm.
// Assumes the length of the array is a power of 2
// Returns the result in complex number space.
func DFT2Radix1D(signals []complex128) ([]complex128, error) {
	if len(signals)  == 0 {
		return make([]complex128, 0), errors.New("DFT2Radix1D: Input array must have size at least one")
	}
	if bits.OnesCount32(uint32(len(signals))) != 1 {
		return make([]complex128, 0), errors.New("DFT2Radix1D: Input array must have size a power of two")
	}
	length := uint32(len(signals))
	numBits := uint32(32-(bits.LeadingZeros32(length)+1))
	shift := 32 - numBits

	// Bit reversal
	for i := uint32(0); i < length; i++ {
		j := bits.Reverse32(i) >> shift
		if j > i {
			signals[i], signals[j] = signals[j], signals[i]
		}
	}

	// radix-2 butterfly
	for window := uint32(2); window <= length; window *= 2 {
		halfWindow := window / 2
		for start := uint32(0); start < length; start += window {
			k := uint32(0)
			for first := start; first < start + halfWindow; first++ {
				second := first + halfWindow
				w := complex(math.Cos(float64(-2) * float64(k)* math.Pi / float64(window)),
					math.Sin(float64(-2) * float64(k)* math.Pi / float64(window)))
				term := w * signals[second]
				signals[second] = signals[first] - term
				signals[first] = signals[first] + term
				k++
			}
		}
	}

	return signals, nil
}

// DFTInverse2Radix1D computes the inverse discrete fourier transform of the given array in the complex number space.
// The calculation is done in place using Cooley-Tukey radix-2 algorithm.
// Assumes the length of the array is a power of 2
// Returns the result in complex number space.
func DFTInverse2Radix1D(signals []complex128) ([]complex128, error) {
	if len(signals)  == 0 {
		return make([]complex128, 0), errors.New("DFT2Radix1D: Input array must have size at least one")
	}
	if bits.OnesCount32(uint32(len(signals))) != 1 {
		return make([]complex128, 0), errors.New("DFT2Radix1D: Input array must have size a power of two")
	}
	inverseSignals := make([]complex128, len(signals))
	for i := 0; i < len(signals); i++ {
		inverseSignals[i] = complex(imag(signals[i]),real(signals[i]))
	}
	result, err := DFT2Radix1D(inverseSignals)
	N := float64(len(signals))
	for i, signal := range result {
		newSignal := complex(imag(signal)/N, real(signal)/N)
		result[i] = newSignal
	}
	return result, err
}

// DFT2Radix1DReal computes the discrete fourier transform of the given array in the real number space.
// The calculation is done in place using Cooley-Tukey radix-2 algorithm.
// Assumes the length of the array is a power of 2
// Returns the result in complex number space.
func DFT2Radix1DReal(signals []float64) ([]complex128, error) {
	complexSignals := make([]complex128, len(signals))
	for i, signal := range signals {
		complexSignals[i] = complex(signal, 0.0)
	}
	return DFT2Radix1D(complexSignals)
}

// DFTInverse2Radix1DReal computes the inverse discrete fourier transform of the given array in the complex number space.
// The calculation is done in place using Cooley-Tukey radix-2 algorithm.
// Assumes the length of the array is a power of 2
// Returns the result in real number space.
func DFTInverse2Radix1DReal(signals []complex128) ([]float64, error) {
	complexSignals, err := DFTInverse2Radix1D(signals)
	realSignals := make([]float64, len(signals))
	for i, signal := range complexSignals {
		realSignals[i] = real(signal)
	}
	return realSignals, err
}

// DFTNaive2DReal computes the discrete fourier transform of the given 2d-array in the real number space.
// Returns the result in complex number space.
func DFT2Radix2DReal(signals [][]float64) ([][]complex128, error) {
	complexSignals := make([][]complex128, len(signals))
	for i, signal := range signals {
		complexSignals[i] = make([]complex128, len(signal))
		for j, num := range signal {
			complexSignals[i][j] = complex(num,0.0)
		}
	}
	return DFT2Radix2D(complexSignals)
}

// DFTInverseNaive2DReal computes the inverse discrete fourier transform of the given 2d-array in the complex number space.
// Returns the result in real number space.
func DFTInverse2Radix2DReal(signals [][]complex128) ([][]float64, error) {
	complexSignals, err := DFTInverseNaive2D(signals)
	realSignals := make([][]float64, len(signals))
	for i, signal := range complexSignals {
		realSignals[i] = make([]float64, len(signal))
		for j, num := range signal {
			realSignals[i][j] = real(num)
		}
	}
	return realSignals, err
}

// DFTNaive2D computes the discrete fourier transform of the given 2d-array in the complex number space.
// Returns the result in complex number space.
func DFT2Radix2D(signals [][]complex128) ([][]complex128, error) {
	return dft2D(signals, true, "radix2")
}

// DFTInverseNaive2D computes the inverse discrete fourier transform of the given 2d-array in the complex number space.
// Returns the result in complex number space.
func DFTInverseRadix2D(signals [][]complex128) ([][]complex128, error) {
	return dft2D(signals, false, "radix2")
}

// DFTNaive1D computes the discrete fourier transform of the given array in the complex number space.
// Returns the result in complex number space.
func DFTNaive1D(signals []complex128) ([]complex128, error) {
	if len(signals)  == 0 {
		return make([]complex128, 0), errors.New("DFTNaive1D: Input array must have size at least one")
	}
	result := make([]complex128, len(signals))
	for n := 0; n < len(signals); n++ {
		sum := complex(0, 0)
		for k := 0; k < len(signals); k++ {
			w := 2 * math.Pi * float64(k) * float64(n) / float64(len(signals))
			sum += signals[k] * complex(math.Cos(w), -math.Sin(w))
		}
		result[n] = sum
	}
	return result, nil
}

// DFTInverseNaive1D computes the inverse discrete fourier transform of the given array in the complex number space.
// Returns the result in complex number space.
func DFTInverseNaive1D(signals []complex128) ([]complex128, error) {
	if len(signals)  == 0 {
		return make([]complex128, 0), errors.New("DFTInverseNaive1D: Input array must have size at least one")
	}
	result := make([]complex128, len(signals))
	for n := 0; n < len(signals); n++ {
		sum := complex(0, 0)
		for k := 0; k < len(signals); k++ {
			w := 2 * math.Pi * float64(k) * float64(n) / float64(len(signals))
			sum += signals[k] * complex(math.Cos(w), math.Sin(w))
		}
		result[n] = sum / complex(float64(len(signals)), 0)
	}
	return result, nil
}

// DFTNaive1DReal computes the discrete fourier transform of the given array in the real number space.
// Returns the result in complex number space.
func DFTNaive1DReal(signals []float64) ([]complex128, error) {
	complexSignals := make([]complex128, len(signals))
	for i, signal := range signals {
		complexSignals[i] = complex(signal, 0.0)
	}
	return DFTNaive1D(complexSignals)
}

// DFTInverseNaive1DReal computes the inverse discrete fourier transform of the given array in the complex number space.
// Returns the result in real number space.
func DFTInverseNaive1DReal(signals []complex128) ([]float64, error) {
	complexSignals, err := DFTInverseNaive1D(signals)
	realSignals := make([]float64, len(signals))
	for i, signal := range complexSignals {
		realSignals[i] = real(signal)
	}
	return realSignals, err
}

// DFTNaive2DReal computes the discrete fourier transform of the given 2d-array in the real number space.
// Returns the result in complex number space.
func DFTNaive2DReal(signals [][]float64) ([][]complex128, error) {
	complexSignals := make([][]complex128, len(signals))
	for i, signal := range signals {
		complexSignals[i] = make([]complex128, len(signal))
		for j, num := range signal {
			complexSignals[i][j] = complex(num,0.0)
		}
	}
	return DFTNaive2D(complexSignals)
}

// DFTInverseNaive2DReal computes the inverse discrete fourier transform of the given 2d-array in the complex number space.
// Returns the result in real number space.
func DFTInverseNaive2DReal(signals [][]complex128) ([][]float64, error) {
	complexSignals, err := DFTInverseNaive2D(signals)
	realSignals := make([][]float64, len(signals))
	for i, signal := range complexSignals {
		realSignals[i] = make([]float64, len(signal))
		for j, num := range signal {
			realSignals[i][j] = real(num)
		}
	}
	return realSignals, err
}

// DFTNaive2D computes the discrete fourier transform of the given 2d-array in the complex number space.
// Returns the result in complex number space.
func DFTNaive2D(signals [][]complex128) ([][]complex128, error) {
	return dft2D(signals, true, "naive")
}

// DFTInverseNaive2D computes the inverse discrete fourier transform of the given 2d-array in the complex number space.
// Returns the result in complex number space.
func DFTInverseNaive2D(signals [][]complex128) ([][]complex128, error) {
	return dft2D(signals, false, "naive")
}

func dft2D(signals [][]complex128, forward bool, algorithm string) ([][]complex128, error) {
	var err error
	height := len(signals)
	// check that input has at least one row
	if height == 0 {
		return make([][]complex128, 0), errors.New("DFTInverseNaive2D: Input 2d-array must have at least one row")
	}

	width := len(signals[0])
	// check that input has at least one column
	if width == 0 {
		return make([][]complex128, 0), errors.New("DFTInverseNaive2D: Input 2d-array must have at least one column")
	}

	// create the result array
	result := make([][]complex128, len(signals))
	for i := 0; i < height; i++ {
		result[i] = make([]complex128, width)
	}

	// Apply DFT on rows as 1d arrays
	channels := make([]chan bool, height)
	for i := 0; i < height; i++ {
		channels[i] = make(chan bool)
		go func(i int) {
			if forward {
				switch algorithm {
				case "radix2":
					result[i], _ = DFT2Radix1D(signals[i])
				default:
					result[i], _ = DFTNaive1D(signals[i])
				}
			} else {
				switch algorithm {
				case "radix2":
					result[i], _ = DFTInverseNaive1D(signals[i])
				default:
					result[i], _ = DFTInverseNaive1D(signals[i])
				}
				for j := 0; j < width; j++ {
					result[i][j] = result[i][j] / complex(float64(width), 0)
				}
			}

			channels[i] <- true
		}(i)
	}
	// Wait on channels
	for i := 0; i < height; i++ {
		<-channels[i]
	}

	// Transpose the array
	result = transposeComplex(result)

	// Apply DFT on columns as 1d arrays
	channels = make([]chan bool, width)
	for i := 0; i < width; i++ {
		channels[i] = make(chan bool)
		go func(i int) {
			if forward {
				switch algorithm {
				case "radix2":
					result[i], _ = DFT2Radix1D(signals[i])
				default:
					result[i], _ = DFTNaive1D(result[i])
				}
			} else {
				switch algorithm {
				case "radix2":
					result[i], _ = DFTInverseNaive1D(signals[i])
				default:
					result[i], _ = DFTInverseNaive1D(result[i])
				}
				for j := 0; j < height; j++ {
					result[i][j] = result[i][j] / complex(float64(height), 0)
				}
			}
			channels[i] <- true
		}(i)
	}
	// Wait on channels
	for i := 0; i < width; i++ {
		<-channels[i]
	}

	result = transposeComplex(result)
	if err != nil {
		return result, err
	}

	return result, nil
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
