package go_fourier

import (
	"math"
	"math/cmplx"
)



func DCT1D(signals []float64) ([]float64, error) {
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

func DCTInverse1D(signals []float64) ([]float64, error) {
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

//
//func DCT2D(signals [][]float64) ([][]float64, error) {
//	return make([][]float64, len(signals)), nil
//}
//
//func DCTInverse2D(signals [][]float64) ([][]float64, error) {
//	return make([][]float64, len(signals)), nil
//}
