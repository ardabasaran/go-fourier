package go_fourier

import (
	"errors"
	"math"
	"math/bits"
)

// DFT2Radix1D computes the discrete fourier transform of the given array in the complex number space.
// The calculation is done in place using Cooley-Tukey radix-2 algorithm.
// The result is stored in the given array.
// Assumes the length of the array is a power of 2
// Calculates the result in complex number space.
func DFT2Radix1D(signals []complex128) error {
	if len(signals) == 0 {
		return errors.New("DFT2Radix1D: Input array must have size at least one")
	}
	if bits.OnesCount32(uint32(len(signals))) != 1 {
		return errors.New("DFT2Radix1D: Input array must have size a power of two")
	}
	length := uint32(len(signals))
	numBits := uint32(32 - (bits.LeadingZeros32(length) + 1))
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
			for first := start; first < start+halfWindow; first++ {
				second := first + halfWindow
				w := complex(math.Cos(float64(-2)*float64(k)*math.Pi/float64(window)),
					math.Sin(float64(-2)*float64(k)*math.Pi/float64(window)))
				term := w * signals[second]
				signals[second] = signals[first] - term
				signals[first] = signals[first] + term
				k++
			}
		}
	}

	return nil
}

// DFTInverse2Radix1D computes the inverse discrete fourier transform of the given array in the complex number space.
// The calculation is done in place using Cooley-Tukey radix-2 algorithm.
// The result is stored in the given array.
// Assumes the length of the array is a power of 2
// Calculates the result in complex number space.
func DFTInverse2Radix1D(signals []complex128) error {
	if len(signals) == 0 {
		return errors.New("DFT2Radix1D: Input array must have size at least one")
	}
	if bits.OnesCount32(uint32(len(signals))) != 1 {
		return errors.New("DFT2Radix1D: Input array must have size a power of two")
	}
	for i := 0; i < len(signals); i++ {
		signals[i] = complex(imag(signals[i]), real(signals[i]))
	}
	err := DFT2Radix1D(signals)
	N := float64(len(signals))
	for i, signal := range signals {
		signals[i] = complex(imag(signal)/N, real(signal)/N)
	}
	return err
}

// DFT2Radix2D computes the discrete fourier transform of the given 2d-array in the complex number space.
// The result is stored in the given array.
// Assumes the dimensions of the array is a power of 2
// Calculates the result in complex number space.
func DFT2Radix2D(signals [][]complex128) error {
	return dft2D(signals, true, "radix2")
}

// DFTInverse2Radix2D computes the inverse discrete fourier transform of the given 2d-array in the complex number space.
// The result is stored in the given array.
// Assumes the dimensions of the array is a power of 2
// Calculates the result in complex number space.
func DFTInverse2Radix2D(signals [][]complex128) error {
	transposedSignals := transposeComplex(signals)
	err := dft2D(transposedSignals, false, "radix2")
	result := transposeComplex(transposedSignals)
	for i := 0; i < len(result); i++ {
		for j := 0; j < len(result[0]); j++ {
			signals[i][j] = result[i][j]
		}
	}
	return err
}

// DFTNaive1D computes the discrete fourier transform of the given array in the complex number space.
// The result is stored in the given array.
// Calculates the result in complex number space.
func DFTNaive1D(signals []complex128) error {
	if len(signals) == 0 {
		return errors.New("DFTNaive1D: Input array must have size at least one")
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
	for i := 0; i < len(signals); i++ {
		signals[i] = result[i]
	}
	return nil
}

// DFTInverseNaive1D computes the inverse discrete fourier transform of the given array in the complex number space.
// The result is stored in the given array.
// Calculates the result in complex number space.
func DFTInverseNaive1D(signals []complex128) error {
	if len(signals) == 0 {
		return errors.New("DFTInverseNaive1D: Input array must have size at least one")
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
	for i := 0; i < len(signals); i++ {
		signals[i] = result[i]
	}
	return nil
}

// DFTNaive2D computes the discrete fourier transform of the given 2d-array in the complex number space.
// The result is stored in the given array.
// Calculates the result in complex number space.
func DFTNaive2D(signals [][]complex128) error {
	return dft2D(signals, true, "naive")
}

// DFTInverseNaive2D computes the inverse discrete fourier transform of the given 2d-array in the complex number space.
// The result is stored in the given array.
// Calculates the result in complex number space.
func DFTInverseNaive2D(signals [][]complex128) error {
	return dft2D(signals , false, "naive")
}

// DFT2Radix1DReal computes the discrete fourier transform of the given array in the real number space.
// Assumes the length of the array is a power of 2
// Returns the result in complex number space.
func DFT2Radix1DReal(signals []float64) ([]complex128, error) {
	complexSignals := make([]complex128, len(signals))
	for i, signal := range signals {
		complexSignals[i] = complex(signal, 0.0)
	}
	err := DFT2Radix1D(complexSignals)
	return complexSignals, err
}

// DFTInverse2Radix1DReal computes the inverse discrete fourier transform of the given array in the complex number space.
// Assumes the length of the array is a power of 2
// Returns the result in real number space.
func DFTInverse2Radix1DReal(signals []complex128) ([]float64, error) {
	copySignals := make([]complex128, len(signals))
	copy(copySignals, signals)
	err := DFTInverse2Radix1D(copySignals)
	realSignals := make([]float64, len(signals))
	for i, signal := range copySignals {
		realSignals[i] = real(signal)
	}
	return realSignals, err
}

// DFTNaive2DReal computes the discrete fourier transform of the given 2d-array in the real number space.
// Assumes the dimensions of the array is a power of 2
// Returns the result in complex number space.
func DFT2Radix2DReal(signals [][]float64) ([][]complex128, error) {
	complexSignals := make([][]complex128, len(signals))
	for i, signal := range signals {
		complexSignals[i] = make([]complex128, len(signal))
		for j, num := range signal {
			complexSignals[i][j] = complex(num, 0.0)
		}
	}
	err := DFT2Radix2D(complexSignals)
	return complexSignals, err
}

// DFTInverse2Radix2DReal computes the inverse discrete fourier transform of the given 2d-array in the complex number space.
// Assumes the dimensions of the array is a power of 2
// Returns the result in real number space.
func DFTInverse2Radix2DReal(signals [][]complex128) ([][]float64, error) {
	copySignals := make([][]complex128, len(signals))
	for i := 0; i < len(signals); i++ {
		copySignals[i] = make([]complex128, len(signals[0]))
		for j := 0; j < len(signals[0]); j++ {
			copySignals[i][j] = signals[i][j]
		}
	}
	err := DFTInverse2Radix2D(copySignals)
	realSignals := make([][]float64, len(signals))
	for i, signal := range copySignals {
		realSignals[i] = make([]float64, len(signal))
		for j, num := range signal {
			realSignals[i][j] = real(num)
		}
	}
	return realSignals, err
}

// DFTNaive1DReal computes the discrete fourier transform of the given array in the real number space.
// Returns the result in complex number space.
func DFTNaive1DReal(signals []float64) ([]complex128, error) {
	complexSignals := make([]complex128, len(signals))
	for i, signal := range signals {
		complexSignals[i] = complex(signal, 0.0)
	}
	err := DFTNaive1D(complexSignals)
	return complexSignals, err
}

// DFTInverseNaive1DReal computes the inverse discrete fourier transform of the given array in the complex number space.
// Returns the result in real number space.
func DFTInverseNaive1DReal(signals []complex128) ([]float64, error) {
	copySignals := make([]complex128, len(signals))
	copy(copySignals, signals)
	err := DFTInverseNaive1D(copySignals)
	realSignals := make([]float64, len(signals))
	for i, signal := range copySignals {
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
			complexSignals[i][j] = complex(num, 0.0)
		}
	}
	err := DFTNaive2D(complexSignals)
	return complexSignals, err
}

// DFTInverseNaive2DReal computes the inverse discrete fourier transform of the given 2d-array in the complex number space.
// Returns the result in real number space.
func DFTInverseNaive2DReal(signals [][]complex128) ([][]float64, error) {
	copySignals := make([][]complex128, len(signals))
	for i := 0; i < len(signals); i++ {
		copySignals[i] = make([]complex128, len(signals[0]))
		for j := 0; j < len(signals[0]); j++ {
			copySignals[i][j] = signals[i][j]
		}
	}
	err := DFTInverseNaive2D(copySignals)
	realSignals := make([][]float64, len(signals))
	for i, signal := range copySignals {
		realSignals[i] = make([]float64, len(signal))
		for j, num := range signal {
			realSignals[i][j] = real(num)
		}
	}
	return realSignals, err
}
