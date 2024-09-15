package mimic

import (
	"math"
	"math/rand/v2"
	"time"

	"gonum.org/v1/gonum/stat/distuv"
)

type Point struct {
	X, Y      float64
	timestamp time.Duration
}

func (p Point) Timestamp() time.Duration {
	return p.timestamp
}

type Config struct {
	Points    []Point
	Duration  time.Duration
	Noise     float64
	Frequency int
}

// SigmaLognormal calculates the velocity at time t using the Sigma-Lognormal Model (D = 1)
func SigmaLognormal(t, t0, mu, sigma float64) float64 {
	exponent := -(math.Pow(math.Log(t-t0)-mu, 2) / (2 * math.Pow(sigma, 2)))
	return 1 / ((t - t0) * sigma * math.Sqrt(2*math.Pi)) * math.Exp(exponent)
}

// Position of the mouse in 2D space for one stroke
func MousePosition(t, t0, mu, sigma, thetaStart, thetaEnd, D float64) (float64, float64) {
	erf := math.Erf((math.Log(t-t0) - mu) / (sigma * math.Sqrt(2)))
	phi := thetaStart + ((thetaEnd-thetaStart)/2)*(1+erf)
	scale := D / (thetaEnd - thetaStart)
	x := scale * (math.Sin(phi) - math.Sin(thetaStart))
	y := scale * (-math.Cos(phi) + math.Cos(thetaStart))
	return x, y
}

// Calculate the end angle θ_ej (Equation 6)
func CalculateEndAngle(xf, yf, t, t0, mu, sigma, thetaStart float64) float64 {
	erf := math.Erf((math.Sqrt(2) * (mu - math.Log(t-t0))) / (2 * sigma))
	atan := math.Atan((xf*math.Tan(thetaStart/2) - yf) / (xf + yf*math.Tan(thetaStart/2)))
	return (thetaStart*erf + thetaStart + 4*atan) / (erf - 1)
}

// CalculateDistanceScalingFactorD implements Equation 8 from the paper
func CalculateDistanceScalingFactor(xf, yf, t, t0, mu, sigma, thetaStart float64) float64 {
	erf := math.Erf((math.Sqrt(2) * (mu - math.Log(t-t0))) / (2 * sigma))
	atan := math.Atan((xf*math.Tan(thetaStart/2) - yf) / (xf + yf*math.Tan(thetaStart/2)))
	return (xf * (thetaStart*(erf-1) - thetaStart*erf - thetaStart - 4*atan)) / ((math.Sin(thetaStart) + math.Sin(2*atan)) * (erf - 1))
}

// Generate a mouse path using the Sigma-Lognormal model
func Generate(config Config) []Point {
	var points []Point

	t0 := 0.0
	mu := 0.0
	sigma := rand.Float64()*2.9 + 0.1

	var totalDelta float64
	deltas := make([]float64, len(config.Points)-1)
	for i := range deltas {
		deltaX := config.Points[i+1].X - config.Points[i].X
		deltaY := config.Points[i+1].Y - config.Points[i].Y
		deltas[i] = math.Sqrt(deltaX*deltaX + deltaY*deltaY)
		totalDelta += deltas[i]
	}
	durations := make([]float64, len(config.Points)-1)
	for i := range deltas {
		durations[i] = deltas[i] / totalDelta * config.Duration.Seconds()
	}

	for i := 0; i < len(config.Points)-1; i++ {
		x0, y0 := config.Points[i].X, config.Points[i].Y
		xf, yf := config.Points[i+1].X, config.Points[i+1].Y
		tj := durations[i]
		frequency := float64(config.Frequency)

		var thetaStart float64
		switch {
		case xf > x0 && yf > y0:
			thetaStart = rand.Float64() * math.Pi / 2
		case xf < x0 && yf > y0:
			thetaStart = math.Pi/2 + rand.Float64()*math.Pi/2
		case xf < x0 && yf < y0:
			thetaStart = -math.Pi/2 + rand.Float64()*math.Pi/2
		case xf > x0 && yf < y0:
			thetaStart = -rand.Float64() * math.Pi / 2
		}

		D := CalculateDistanceScalingFactor(xf-x0, yf-y0, tj, t0, mu, sigma, thetaStart)
		thetaEnd := CalculateEndAngle(xf-x0, yf-y0, tj, t0, mu, sigma, thetaStart)

		for t := 0.0; t < tj; t += 1 / frequency {
			x, y := MousePosition(t, t0, mu, sigma, thetaStart, thetaEnd, D)
			noise := distuv.Normal{Mu: 0, Sigma: config.Noise}
			x += noise.Rand()
			y += noise.Rand()
			points = append(points, Point{X: x + x0, Y: y + y0, timestamp: time.Duration(t) * time.Second})
		}
	}

	return points
}
