package mimic

import (
	"math"
	"math/rand/v2"
	"time"

	"gonum.org/v1/gonum/stat/distuv"
)

type Point struct {
	X, Y float64

	velocity  float64
	timestamp time.Duration
}

func (p Point) Velocity() float64 {
	return p.velocity
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

func SigmaLognormalDerivative(t, t0, mu, sigma float64) float64 {
	exponent := -(math.Pow(math.Log(t-t0)-mu, 2) / (2 * math.Pow(sigma, 2)))
	common := math.Exp(exponent) / ((t - t0) * sigma * math.Sqrt(2*math.Pi))

	term1 := -1 / (t - t0)
	term2 := (math.Log(t-t0) - mu) / (sigma * sigma * (t - t0))

	return (term1 + term2) * common
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

// Function to constrain the velocity within the thresholds
func ConstrainVelocity(mu, sigma, t0, D float64) (float64, float64) {
	vthreshold := 1500.0
	vthresholdmin := 300.0

	tmax := t0 + math.Exp(-sigma*sigma+mu)
	vmax := D * SigmaLognormal(tmax, t0, mu, sigma)

	for vmax > vthreshold || vmax < vthresholdmin {
		switch {
		case vmax > vthreshold:
			if rand.Float64() < 0.5 {
				mu = AdjustMu(mu, sigma, tmax, t0, vthreshold, D)
			} else {
				D = AdjustD(D, vmax, vthreshold, vthresholdmin)
			}
		case vmax < vthresholdmin:
			D = AdjustD(D, vmax, vthreshold, vthresholdmin)
		}

		tmax = t0 + math.Exp(-sigma*sigma+mu)
		vmax = D * SigmaLognormal(tmax, t0, mu, sigma)
	}

	return mu, D
}

// Function to adjust mu to reduce/increase velocity to fit within the range
func AdjustMu(mu, sigma, t, t0, vmax, D float64) float64 {
	mu1 := -sigma*math.Sqrt(2*math.Log(-(D/(sigma*vmax*(t0-t))))-math.Log(2*math.Pi)) + math.Log(-t0+t)
	mu2 := sigma*math.Sqrt(2*math.Log(D/(sigma*vmax*(-t0+t)))-math.Log(2*math.Pi)) + math.Log(-t0+t)

	derivative1 := SigmaLognormalDerivative(t, t0, mu1, sigma)
	derivative2 := SigmaLognormalDerivative(t, t0, mu2, sigma)

	switch {
	case derivative1 < 0:
		return mu1
	case derivative2 < 0:
		return mu2
	default:
		return mu
	}
}

// Function to adjust D (distance scaling factor) to fit the velocity within the range
func AdjustD(D, vmax, vthreshold, vthresholdmin float64) float64 {
	dmin := (vthresholdmin * D) / vmax
	dmax := (vthreshold * D) / vmax

	return rand.Float64()*(dmax-dmin) + dmin
}

// Generate a mouse path using the Sigma-Lognormal model
func Generate(config Config) []Point {
	var points []Point

	t0 := 0.0
	frequency := float64(config.Frequency)

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

		var mu, sigma, D, thetaStart, thetaEnd float64
		for {
			thetaStart = math.Atan2(yf-y0, xf-x0)
			thetaStart += (rand.Float64() - 0.5) * (math.Pi / 4)
			mu = rand.Float64()*6 - 3
			sigma = rand.Float64()*2.9 + 0.1
			D = rand.Float64()*1100 + 100

			mu, D = ConstrainVelocity(mu, sigma, t0, D)
			D = CalculateDistanceScalingFactor(xf-x0, yf-y0, tj, t0, mu, sigma, thetaStart)
			thetaEnd = CalculateEndAngle(xf-x0, yf-y0, tj, t0, mu, sigma, thetaStart)
			if !(math.IsNaN(D) || math.IsNaN(thetaEnd) || math.IsInf(D, 0) || math.IsInf(thetaEnd, 0)) {
				break
			}
		}

		for t := 0.0; t < tj; t += 1 / frequency {
			x, y := MousePosition(t, t0, mu, sigma, thetaStart, thetaEnd, D)
			noise := distuv.Normal{Mu: 0, Sigma: config.Noise}
			x += noise.Rand()
			y += noise.Rand()
			v := SigmaLognormal(t, t0, mu, sigma)
			points = append(points, Point{X: x + x0, Y: y + y0, velocity: v, timestamp: time.Duration(t) * time.Second})
		}
	}

	return points
}
