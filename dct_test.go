package go_fourier

import (
	"testing"
)

var testsDCT = []struct {
	input  []float64
	output []float64
}{
	{
		[]float64{4., 7., 2., 5., 6., 9., 1., 3.},
		[]float64{13.08147545, 0.9427605, -2.42178421, 3.54099744, -0.35355339, -3.76312166, 0.62045243, -3.98891654},
	},
	{
		[]float64{0.4, 0.2, 0.75, 0.11},
		[]float64{0.73, 0.0406227, -0.22, 0.43777825},
	},
}

var testsDCT2D = [] struct{
	input [][]float64
	output [][]float64
}{
	{
		[][]float64{{0.1, 0.4, 2.5, 15}, {1.35, 20.234, 6.324, 0.532}, {4.12, 3.12, 2.12, 1.12}, {9, 8, 7, 6}},
		[][]float64{{21.73, -0.77143093, -3.119, -5.60439536}, {-1.4897184, -7.62772449, 0.6463783, -4.28325548}, {2.27, -7.30027656, 9.219, 3.10288819}, {-7.49005601, -4.74225548, 9.71083504, 5.13272449}},
	},
	{
		[][]float64{{1, 2}, {3, 4}},
		[][]float64{{5, -1}, {-2, 0}},
	},
	{
		[][]float64{{1, 2}, {3, 4}, {5, 6}, {7, 8}},
		[][]float64{{12.72792206, -1.41421356}, {-6.30864406, 0}, {0, 0}, {-0.44834153, 0}},
	},
}

func TestDCT1D(t *testing.T) {
	for _, test := range testsDCT {
		actual, _ := DCT1D(test.input)
		expected := test.output
		for i, c := range expected {
			diff := actual[i] - c
			if diff > 1e-2 {
				t.Errorf("Difference of %v and %v is %v", actual[i], c, diff)
			}
		}
	}
}

func TestDCTInverse1D(t *testing.T) {
	for _, test := range testsDCT {
		actual, _ := DCTInverse1D(test.output)
		expected := test.input
		for i, c := range expected {
			diff := actual[i] - c
			if diff > 1e-2 {
				t.Errorf("Difference of %v and %v is %v", actual[i], c, diff)
			}
		}
	}
}

func TestDCT2D(t *testing.T) {
	for _, test := range testsDCT2D {
		actual, _ := DCT2D(test.input)
		expected := test.output
		for i := 0; i < len(test.input); i++ {
			for j := 0; j < len(test.input[0]); j++ {
				diff := actual[i][j] - expected[i][j]
				if diff > 1e-2 {
					t.Errorf("Difference of %v and %v is %v", actual[i][j], expected[i][j], diff)
				}
			}
		}
	}
}

func TestDCTInverse2D(t *testing.T) {
	for _, test := range testsDCT2D {
		actual, _ := DCTInverse2D(test.output)
		expected := test.input
		for i := 0; i < len(test.output); i++ {
			for j := 0; j < len(test.output[0]); j++ {
				diff := actual[i][j] - expected[i][j]
				if diff > 1e-2 {
					t.Errorf("Difference of %v and %v is %v", actual[i][j], expected[i][j], diff)
				}
			}
		}
	}
}

