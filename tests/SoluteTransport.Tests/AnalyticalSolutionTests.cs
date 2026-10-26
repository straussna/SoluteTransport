using SoluteTransport;

namespace SoluteTransport.Tests;

/// <summary>
/// Compares <see cref="Stream"/> with closed-form solutions of the OTIS
/// equations (Runkel 1998; Runkel 1996) on reaches long enough that the
/// downstream boundary does not affect the compared segments.
/// </summary>
public class AnalyticalSolutionTests
{
    private const double A = 1.0;
    private const double Q = 0.5;
    private const double U = Q / A;
    private const double D = 0.5;
    private const double Dx = 0.5;
    private const double Dt = 0.5;

    private static StreamParameters Conservative() =>
        new(A: A, AS: 0.1, CL: 0.0, D: D, Q: Q, QLin: 0.0,
            Alpha: 0.0, ChatS: 0.0, Kd: 0.0, Lambda: 0.0, LambdaS: 0.0,
            LambdaHat: 0.0, LambdaHatS: 0.0, Rho: 0.0);

    private static double X(int segment) => (segment + 0.5) * Dx;

    private static double[] Pulse(double concentration, double duration, double totalDuration)
    {
        int levels = (int)Math.Round(totalDuration / Dt) + 1;
        return Enumerable.Range(0, levels)
            .Select(j => j * Dt <= duration ? concentration : 0.0)
            .ToArray();
    }

    /// <summary>
    /// exp(a) · erfc(z), evaluated without overflow (Numerical Recipes erfcc,
    /// fractional error below 1.2e-7).
    /// </summary>
    private static double ExpErfc(double a, double z)
    {
        double t = 1 / (1 + 0.5 * Math.Abs(z));
        double poly = -1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418
            + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587
            + t * (-0.82215223 + t * 0.17087277))))))));
        double scaled = t * Math.Exp(a - z * z + poly);
        return z >= 0 ? scaled : 2 * Math.Exp(a) - scaled;
    }

    /// <summary>
    /// Continuous load of concentration <paramref name="c0"/> and duration
    /// <paramref name="tau"/> with first-order decay <paramref name="lambda"/>,
    /// Runkel (1996) equations (4) and (11).
    /// </summary>
    private static double FiniteDurationLoad(double x, double t, double c0, double tau, double lambda)
    {
        double gammaR = Math.Sqrt(1 + 2 * (2 * lambda * D / (U * U)));
        double minus = U * x / (2 * D) * (1 - gammaR);
        double plus = U * x / (2 * D) * (1 + gammaR);

        double Term(double time) => time <= 0
            ? 0
            : ExpErfc(minus, (x - U * time * gammaR) / (2 * Math.Sqrt(D * time)))
              + ExpErfc(plus, (x + U * time * gammaR) / (2 * Math.Sqrt(D * time)));

        return c0 / 2 * (Term(t) - Term(t - tau));
    }

    /// <summary>
    /// Steady state of -U C' + D C'' + k (Cinf - C) = 0 with C(0) = Cbc on a
    /// semi-infinite reach.
    /// </summary>
    private static double SteadyState(double x, double cbc, double k, double cInf)
    {
        double r = (U - Math.Sqrt(U * U + 4 * k * D)) / (2 * D);
        return cInf + (cbc - cInf) * Math.Exp(r * x);
    }

    [Theory]
    [InlineData(0.0, 250.0)]   // conservative, load still on
    [InlineData(0.0, 60.0)]    // conservative, load switched off
    [InlineData(1e-3, 60.0)]   // with first-order decay
    public void FiniteDurationLoad_MatchesRunkel1996(double lambda, double tau)
    {
        const double length = 300, duration = 200, c0 = 1.0;
        var p = Conservative() with { Lambda = lambda };
        var result = new Stream(length, duration, Dt, Dx, Pulse(c0, tau, duration), p)
            .SimulateAdvectionAndDispersion();

        int last = result.C.GetLength(1) - 1;
        double t = last * Dt;
        double maxError = 0;
        for (int i = 0; X(i) <= 150; i++)
        {
            maxError = Math.Max(maxError,
                Math.Abs(result.C[i, last] - FiniteDurationLoad(X(i), t, c0, tau, lambda)));
        }

        Assert.True(maxError < 0.01, $"max |C - analytical| = {maxError}");
    }

    [Fact]
    public void ZeroUpstreamConcentration_StaysZero()
    {
        var p = Conservative() with { Alpha = 1e-3, Kd = 0.1, LambdaHat = 1e-3, Rho = 1.0 };
        var result = new Stream(50, 100, Dt, Dx, [], p).SimulateAdvectionAndDispersion();

        foreach (double value in result.C) Assert.Equal(0.0, value);
        foreach (double value in result.CS) Assert.Equal(0.0, value);
        foreach (double value in result.CSed) Assert.Equal(0.0, value);
    }

    [Fact]
    public void ConstantInput_WithStorageAndSorption_ReachesUniformEquilibrium()
    {
        const double cbc = 2.0, kd = 0.3;
        var p = Conservative() with { Alpha = 0.01, AS = 0.5, Kd = kd, LambdaHat = 0.01, Rho = 2.0 };
        var result = new Stream(50, 10000, Dt, Dx, Pulse(cbc, double.MaxValue, 10000), p)
            .SimulateAdvectionAndDispersion();

        int last = result.C.GetLength(1) - 1;
        for (int i = 0; i < result.C.GetLength(0); i++)
        {
            Assert.Equal(cbc, result.C[i, last], 6);
            Assert.Equal(cbc, result.CS[i, last], 6);
            Assert.Equal(kd * cbc, result.CSed[i, last], 6);
        }
    }

    [Theory]
    [InlineData("decay")]
    [InlineData("storage decay")]
    [InlineData("lateral inflow")]
    public void ConstantInput_ReachesSteadyStateProfile(string process)
    {
        const double cbc = 1.0, alpha = 0.01, aS = 0.5, lambdaS = 0.005;
        var p = Conservative();
        double k, cInf = 0;
        switch (process)
        {
            case "decay":
                p = p with { Lambda = 0.002 };
                k = 0.002;
                break;
            case "storage decay":
                // Runkel (1998) equation (32): effective main-channel decay.
                p = p with { Alpha = alpha, AS = aS, LambdaS = lambdaS };
                k = alpha * aS * lambdaS / (alpha * A + aS * lambdaS);
                break;
            default:
                p = p with { QLin = 0.001, CL = 0.4 };
                k = 0.001 / A;
                cInf = 0.4;
                break;
        }

        var result = new Stream(400, 6000, Dt, Dx, Pulse(cbc, double.MaxValue, 6000), p)
            .SimulateAdvectionAndDispersion();

        int last = result.C.GetLength(1) - 1;
        double maxError = 0;
        for (int i = 0; X(i) <= 200; i++)
        {
            maxError = Math.Max(maxError, Math.Abs(result.C[i, last] - SteadyState(X(i), cbc, k, cInf)));
        }

        Assert.True(maxError < 1e-3, $"max |C - steady state| = {maxError}");
    }

    [Fact]
    public void DownstreamDispersiveFlux_MatchesSteadyStateProfile()
    {
        const double length = 50, cbc = 1.0, flux = -0.01;
        var result = new Stream(length, 4000, Dt, Dx, Pulse(cbc, double.MaxValue, 4000), Conservative(), flux)
            .SimulateAdvectionAndDispersion();

        // -U C' + D C'' = 0, C(0) = Cbc, D C'(L) = flux.
        double Exact(double x) =>
            cbc + flux / U * (Math.Exp(U * (x - length) / D) - Math.Exp(-U * length / D));

        int last = result.C.GetLength(1) - 1;
        double maxError = 0;
        for (int i = 0; i < result.C.GetLength(0); i++)
        {
            maxError = Math.Max(maxError, Math.Abs(result.C[i, last] - Exact(X(i))));
        }

        Assert.True(maxError < 1e-3, $"max |C - steady state| = {maxError}");
    }

    [Fact]
    public void TimeStepRefinement_ConvergesAtSecondOrder()
    {
        // Smooth upstream input so the Crank-Nicolson error is O(dt²).
        const double length = 50, duration = 400;
        var p = Conservative() with
        {
            Alpha = 0.01, AS = 0.5, Lambda = 1e-3, LambdaS = 2e-3,
            Kd = 0.2, LambdaHat = 0.02, Rho = 1.5,
        };

        double[,] Run(double dt)
        {
            int levels = (int)Math.Round(duration / dt) + 1;
            double[] boundary = Enumerable.Range(0, levels)
                .Select(j => 1 - Math.Exp(-Math.Pow(j * dt / 50, 2)))
                .ToArray();
            return new Stream(length, duration, dt, Dx, boundary, p).SimulateAdvectionAndDispersion().C;
        }

        double[,] coarse = Run(2), medium = Run(1), fine = Run(0.5);
        double e1 = 0, e2 = 0;
        for (int i = 0; i < coarse.GetLength(0); i++)
        {
            e1 = Math.Max(e1, Math.Abs(coarse[i, 200] - medium[i, 400]));
            e2 = Math.Max(e2, Math.Abs(medium[i, 400] - fine[i, 800]));
        }

        double order = Math.Log2(e1 / e2);
        Assert.InRange(order, 1.8, 2.2);
    }
}
