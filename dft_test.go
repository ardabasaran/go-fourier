package go_fourier

import (
	"math"
	"testing"
)

var tests = []struct {
	input  []float64
	output []complex128
}{
	{
		[]float64{11.0, 22.0, 3.0, 4.0},
		[]complex128{complex(40.0, 0.0), complex(8, -18), complex(-12, 0), complex(8, 18)},
	},
	{
		[]float64{1.0, 2.0, 3.0, 4.0, 5.0},
		[]complex128{complex(15.0, 0.0), complex(-2.5, 3.441), complex(-2.5, 0.812),
			complex(-2.5, -0.812), complex(-2.5, -3.441)},
	},
	{
		[]float64{1.5, 10.2, 34.5, 102.1, 40.9, 3.25},
		[]complex128{complex(192.45, 0.0), complex(-131.575, -0.476), complex(59.175, -11.561),
			complex(-38.65, 0), complex(59.175, 11.561), complex(-131.575, 0.476)},
	},
}

func TestDFTNaive1DReal(t *testing.T) {
	for _, test := range tests {
		actual := DFTNaive1DReal(test.input)
		expected := test.output
		t.Log(actual)
		for i, c := range actual {
			t.Logf("Expected: %v \tActual:%v\n", expected[i], c)
			imagDiffSq := (imag(expected[i]) - imag(c)) * (imag(expected[i]) - imag(c))
			realDiffSq := (real(expected[i]) - real(c)) * (real(expected[i]) - real(c))
			sqrtDiff := math.Sqrt(imagDiffSq + realDiffSq)
			if sqrtDiff > 1e-2 {
				t.Errorf("Difference of %v and %v is %v", expected[i], c, sqrtDiff)
			}
		}
	}
}

func TestDFTInverseNaive1DReal(t *testing.T) {
	for _, test := range tests {
		actual := DFTInverseNaive1DReal(test.output)
		expected := test.input
		for i, c := range actual {
			realDiffSq := (expected[i] - c) * (expected[i] - c)
			sqrtDiff := math.Sqrt(realDiffSq)
			if sqrtDiff > 1e-2 {
				t.Errorf("Difference of %v and %v is %v", expected[i], c, sqrtDiff)
			}
		}
	}
}
