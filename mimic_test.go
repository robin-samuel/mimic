package mimic_test

import (
	"mimic"
	"testing"
	"time"

	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
)

func TestMimic(t *testing.T) {
	config := mimic.Config{
		Points: []mimic.Point{
			{X: 400, Y: 100},
			{X: 500, Y: 500},
			{X: 600, Y: 400},
		},
		Duration:  2 * time.Second,
		Noise:     0.2,
		Frequency: 60,
	}
	path := mimic.Generate(config)

	points := make(plotter.XYs, 0)
	for _, p := range path {
		points = append(points, plotter.XY{X: p.X, Y: p.Y})
	}

	p := plot.New()
	p.Title.Text = "Mouse Path"
	p.X.Label.Text = "X"
	p.Y.Label.Text = "Y"

	s, err := plotter.NewScatter(points)
	if err != nil {
		panic(err)
	}
	s.GlyphStyle.Radius = vg.Points(1)
	p.Add(s)

	if err := p.Save(10*vg.Inch, 10*vg.Inch, "mouse_path.png"); err != nil {
		panic(err)
	}
}
