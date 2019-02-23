package go_fourier

import (
	"math"
)

func DFTNaive1D(signals []complex128) []complex128 {
	result := make([]complex128, len(signals))
	for n := 0; n < len(signals); n++ {
		sum := complex(0, 0)
		for k := 0; k < len(signals); k++ {
			w := 2 * math.Pi * float64(k) * float64(n) / float64(len(signals))
			sum += signals[k] * complex(math.Cos(w), -math.Sin(w))
		}
		result[n] = sum
	}
	return result
}

func DFTInverseNaive1D(signals []complex128) []complex128 {
	result := make([]complex128, len(signals))
	for n := 0; n < len(signals); n++ {
		sum := complex(0, 0)
		for k := 0; k < len(signals); k++ {
			w := 2 * math.Pi * float64(k) * float64(n) / float64(len(signals))
			sum += signals[k] * complex(math.Cos(w), math.Sin(w))
		}
		result[n] = sum / complex(float64(len(signals)), 0)
	}
	return result
}

func DFTNaive1DReal(signals []float64) []complex128 {
	complexSignals := make([]complex128, len(signals))
	for i, signal := range signals {
		complexSignals[i] = complex(signal, 0.0)
	}
	return DFTNaive1D(complexSignals)
}

func DFTInverseNaive1DReal(signals []complex128) []float64 {
	complexSignals := DFTInverseNaive1D(signals)
	realSignals := make([]float64, len(signals))
	for i, signal := range complexSignals {
		realSignals[i] = real(signal)
	}
	return realSignals
}
