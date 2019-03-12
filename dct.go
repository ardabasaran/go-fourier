package go_fourier

import (
	"errors"
	"math"
	"math/bits"
	"math/cmplx"
)

// DCT1D computes the discrete cosine transform of the given array in the complex number space.
// Assumes the length of the array is a power of 2
// Returns the result in real number space.
func DCT1D(signals []float64) ([]float64, error) {
	if len(signals) == 0 {
		return make([]float64,0), errors.New("DCT1D: Input array must have size at least one")
	}
	if bits.OnesCount32(uint32(len(signals))) != 1 {
		return make([]float64,0),  errors.New("DCT1D: Input array must have size a power of two")
	}

	N := len(signals)
	y := make([]complex128, N)
	for i := 0; i < N/2; i++ {
		y[i] = complex(signals[2*i],0.0)
		y[N-1-i] = complex(signals[2*i+1], 0.0)
	}
	err := DFT2Radix1D(y)
	result := make([]float64, len(signals))
	sqrtTermForFirst := math.Sqrt(1.0/(float64(N)))
	sqrtTermForRest := math.Sqrt(2.0/(float64(N)))
	for n := 0; n < N; n++ {
		shift := cmplx.Exp(-1i*math.Pi*complex(float64(n)/float64(2*N),0))
		result[n] = real(y[n]*shift)
		if n == 0 {
			result[n] *= sqrtTermForFirst
		} else {
			result[n] *= sqrtTermForRest
		}
	}
	return result, err
}

// DCTInverse1D computes the inverse discrete cosine transform of the given array in the complex number space.
// Assumes the length of the array is a power of 2
// Returns the result in complex number space.
func DCTInverse1D(signals []float64) ([]float64, error) {
	if len(signals) == 0 {
		return make([]float64,0), errors.New("DCTInverse1D: Input array must have size at least one")
	}
	if bits.OnesCount32(uint32(len(signals))) != 1 {
		return make([]float64,0),  errors.New("DCTInverse1D: Input array must have size a power of two")
	}

	N := len(signals)
	complexSignals := make([]complex128, len(signals))
	for n := 0; n < N; n++ {
		shift := cmplx.Exp(1i*math.Pi*complex(float64(n)/float64(2*N),0))
		complexSignals[n] = complex(signals[n] * math.Sqrt(2.0/(float64(N))), 0.0) * shift
	}
	complexSignals[0] /= complex(math.Sqrt(2.0), 0.0)

	err := DFTInverse2Radix1D(complexSignals)
	result := make([]float64, len(signals))
	for i := 0; i < N/2; i++ {
		result[2*i] = float64(N)*real(complexSignals[i])
		result[2*i+1] = float64(N)*real(complexSignals[N-1-i])
	}
	return result, err
}

// DCT2D computes the discrete cosine transform of the given 2d-array in the complex number space.
// Assumes the dimensions of the array is a power of 2
// Returns the result in complex number space.
func DCT2D(signals [][]float64) ([][]float64, error) {
	return dct2D(signals, true)
}

// DCTInverse2D computes the inverse discrete cosine transform of the given 2d-array in the complex number space.
// Assumes the dimensions of the array is a power of 2
// Returns the result in complex number space.
func DCTInverse2D(signals [][]float64) ([][]float64, error) {
	transposedSignals := transposeReal(signals)
	result, err := dct2D(transposedSignals, false)
	result = transposeReal(result)
	return result, err
}


