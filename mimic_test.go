package mimic_test

import (
	"testing"
	"time"

	"github.com/robin-samuel/mimic"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
)

// Conversion factor: pixels per inch
const dpi = 96.0

func TestMimic(t *testing.T) {
	config := mimic.Config{
		Points: []mimic.Point{
			{X: 100, Y: 100},
			{X: 500, Y: 500},
			{X: 600, Y: 400},
			{X: 300, Y: 1000},
			{X: 50, Y: 600},
		},
		Duration:  4 * time.Second,
		Noise:     0.05,
		Frequency: 60,
		Viewport: &mimic.Viewport{
			Width:  1920,
			Height: 1080 - 95,
		},
	}
	path := mimic.Generate(config)

	points := make(plotter.XYs, 0)
	for _, p := range path {
		points = append(points, plotter.XY{X: p.X, Y: p.Y})
	}

	p := plot.New()
	p.X.Min = 0
	p.X.Max = float64(config.Viewport.Width)
	p.Y.Min = 0
	p.Y.Max = float64(config.Viewport.Height)

	s, err := plotter.NewScatter(points)
	if err != nil {
		panic(err)
	}
	s.GlyphStyle.Radius = vg.Points(1)
	p.Add(s)

	width := float64(config.Viewport.Width) / dpi
	height := float64(config.Viewport.Height) / dpi

	if err := p.Save(vg.Length(width)*vg.Inch, vg.Length(height)*vg.Inch, "sample.png"); err != nil {
		panic(err)
	}
}
