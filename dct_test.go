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
		[]float64{13.08147545,  0.9427605 , -2.42178421,  3.54099744, -0.35355339, -3.76312166,  0.62045243, -3.98891654},
	},
	{
		[]float64{0.4, 0.2, 0.75, 0.11},
		[]float64{0.73, 0.0406227, -0.22, 0.43777825},
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

//func TestDCT(t *testing.T)  {
//	for _, test := range testsDCT {
//		fmt.Println(test.input)
//		transformed, _ := DCT1D(test.input)
//		fmt.Println(transformed)
//		inverse, _ := DCTInverse1D(transformed)
//		fmt.Println(inverse)
//	}
//}
