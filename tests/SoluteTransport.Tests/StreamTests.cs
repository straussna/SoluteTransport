using SoluteTransport;

namespace SoluteTransport.Tests;

public class StreamParametersTests
{
    [Fact]
    public void Constructor_SetsAllProperties()
    {
        var p = new StreamParameters(
            A: 1.0, AS: 0.5, CL: 0.1, D: 0.05, Q: 1.0, QLin: 0.0,
            Alpha: 0.001, ChatS: 0.0, Kd: 0.0, Lambda: 0.0, LambdaS: 0.0,
            LambdaHat: 0.0, LambdaHatS: 0.0, Rho: 0.0);

        Assert.Equal(1.0, p.A);
        Assert.Equal(0.5, p.AS);
        Assert.Equal(0.1, p.CL);
        Assert.Equal(0.05, p.D);
        Assert.Equal(1.0, p.Q);
        Assert.Equal(0.0, p.QLin);
        Assert.Equal(0.001, p.Alpha);
    }
}

public class StreamTests
{
    // ------------------------------------------------------------------
    // Helper: builds a minimal valid StreamParameters record.
    // ------------------------------------------------------------------
    private static StreamParameters DefaultParameters() =>
        new(A: 1.0, AS: 0.1, CL: 0.0, D: 0.05, Q: 0.5, QLin: 0.0,
            Alpha: 0.001, ChatS: 0.0, Kd: 0.0, Lambda: 0.0, LambdaS: 0.0,
            LambdaHat: 0.0, LambdaHatS: 0.0, Rho: 0.0);

    // ------------------------------------------------------------------
    // Construction / argument validation
    // ------------------------------------------------------------------

    [Fact]
    public void Constructor_ThrowsArgumentNullException_WhenUpstreamDataIsNull()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new Stream(10, 10, 1, 1, null!, DefaultParameters()));
    }

    [Fact]
    public void Constructor_ThrowsArgumentNullException_WhenParametersIsNull()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new Stream(10, 10, 1, 1, [], null!));
    }

    [Theory]
    [InlineData(0)]
    [InlineData(-1)]
    public void Constructor_ThrowsArgumentException_WhenDxIsNotPositive(double dx)
    {
        Assert.Throws<ArgumentException>(() =>
            new Stream(10, 10, 1, dx, [], DefaultParameters()));
    }

    [Theory]
    [InlineData(0)]
    [InlineData(-1)]
    public void Constructor_ThrowsArgumentException_WhenDtIsNotPositive(double dt)
    {
        Assert.Throws<ArgumentException>(() =>
            new Stream(10, 10, dt, 1, [], DefaultParameters()));
    }

    [Fact]
    public void Constructor_ThrowsArgumentException_WhenAIsZero()
    {
        var p = DefaultParameters() with { A = 0 };
        Assert.Throws<ArgumentException>(() =>
            new Stream(10, 10, 1, 1, [], p));
    }

    [Fact]
    public void Constructor_ThrowsArgumentException_WhenASIsZero()
    {
        var p = DefaultParameters() with { AS = 0 };
        Assert.Throws<ArgumentException>(() =>
            new Stream(10, 10, 1, 1, [], p));
    }

    // ------------------------------------------------------------------
    // SimulateAdvectionAndDispersion – basic sanity checks
    // ------------------------------------------------------------------

    [Fact]
    public void Simulate_ReturnsResult_WithExpectedDimensions()
    {
        double streamLength = 10.0;
        double totalDuration = 10.0;
        double dt = 1.0;
        double dx = 1.0;

        int numSegments = (int)(streamLength / dx);
        int numTimeSteps = (int)(totalDuration / dt);

        double[] boundary = new double[numTimeSteps];
        var stream = new Stream(streamLength, totalDuration, dt, dx, boundary, DefaultParameters());

        SimulationResult result = stream.SimulateAdvectionAndDispersion();

        Assert.Equal(numSegments, result.C.GetLength(0));
        Assert.Equal(numTimeSteps, result.C.GetLength(1));
        Assert.Equal(numSegments, result.CS.GetLength(0));
        Assert.Equal(numTimeSteps, result.CS.GetLength(1));
        Assert.Equal(numSegments, result.CSed.GetLength(0));
        Assert.Equal(numTimeSteps, result.CSed.GetLength(1));
    }

    [Fact]
    public void Simulate_AllZeroConcentrations_WhenNoInputAndNoTransport()
    {
        // With Q=0 and D=0 the advective-diffusive flux coefficients G and E
        // are both zero, so the simulation with all-zero inputs and all-zero
        // initial conditions must remain identically zero.
        var inertParameters = new StreamParameters(
            A: 1.0, AS: 0.1, CL: 0.0, D: 0.0, Q: 0.0, QLin: 0.0,
            Alpha: 0.0, ChatS: 0.0, Kd: 0.0, Lambda: 0.0, LambdaS: 0.0,
            LambdaHat: 0.0, LambdaHatS: 0.0, Rho: 0.0);

        double[] boundary = new double[100]; // all zeros
        var stream = new Stream(10, 100, 1, 1, boundary, inertParameters);

        SimulationResult result = stream.SimulateAdvectionAndDispersion();

        for (int i = 0; i < result.C.GetLength(0); i++)
        {
            for (int j = 0; j < result.C.GetLength(1); j++)
            {
                Assert.Equal(0.0, result.C[i, j]);
                Assert.Equal(0.0, result.CS[i, j]);
                Assert.Equal(0.0, result.CSed[i, j]);
            }
        }
    }

    [Fact]
    public void Simulate_UpstreamBoundaryConcentration_PropagatesDownstream()
    {
        double streamLength = 10.0;
        double totalDuration = 200.0;
        double dt = 0.5;
        double dx = 1.0;

        int numTimeSteps = (int)(totalDuration / dt);

        // Apply a constant upstream concentration of 1 kg/m³ for the full simulation
        double[] boundary = Enumerable.Repeat(1.0, numTimeSteps).ToArray();

        var stream = new Stream(streamLength, totalDuration, dt, dx, boundary, DefaultParameters());
        SimulationResult result = stream.SimulateAdvectionAndDispersion();

        // The concentration at the last segment at the final time step should
        // be positive because the tracer has had time to travel downstream.
        int lastSegment = (int)(streamLength / dx) - 1;
        int lastTime = numTimeSteps - 1;
        Assert.True(result.C[lastSegment, lastTime] > 0.0,
            $"Expected concentration > 0 at downstream end, got {result.C[lastSegment, lastTime]}");
    }

    [Fact]
    public void Simulate_ConcentrationRemainsFinite()
    {
        double streamLength = 20.0;
        double totalDuration = 100.0;
        double dt = 0.5;
        double dx = 1.0;

        int numTimeSteps = (int)(totalDuration / dt);
        double[] boundary = Enumerable.Repeat(1.0, numTimeSteps).ToArray();

        var stream = new Stream(streamLength, totalDuration, dt, dx, boundary, DefaultParameters());
        SimulationResult result = stream.SimulateAdvectionAndDispersion();

        for (int i = 0; i < result.C.GetLength(0); i++)
        {
            for (int j = 0; j < result.C.GetLength(1); j++)
            {
                Assert.False(double.IsNaN(result.C[i, j]), $"C[{i},{j}] is NaN");
                Assert.False(double.IsInfinity(result.C[i, j]), $"C[{i},{j}] is Infinity");
            }
        }
    }

    [Fact]
    public void Simulate_ConcentrationRemainsNonNegative_WithPositiveUpstreamInput()
    {
        double streamLength = 20.0;
        double totalDuration = 100.0;
        double dt = 0.5;
        double dx = 1.0;

        int numTimeSteps = (int)(totalDuration / dt);
        double[] boundary = Enumerable.Repeat(1.0, numTimeSteps).ToArray();

        var stream = new Stream(streamLength, totalDuration, dt, dx, boundary, DefaultParameters());
        SimulationResult result = stream.SimulateAdvectionAndDispersion();

        for (int i = 0; i < result.C.GetLength(0); i++)
        {
            for (int j = 0; j < result.C.GetLength(1); j++)
            {
                Assert.True(result.C[i, j] >= 0.0,
                    $"C[{i},{j}] = {result.C[i, j]} is negative");
            }
        }
    }

    [Fact]
    public void Simulate_ShortUpstreamBoundary_PadsWithZero()
    {
        double streamLength = 10.0;
        double totalDuration = 100.0;
        double dt = 1.0;
        double dx = 1.0;
        int numTimeSteps = (int)(totalDuration / dt);

        // Only supply 10 time steps of boundary data (less than numTimeSteps)
        double[] boundary = Enumerable.Repeat(1.0, 10).ToArray();

        var stream = new Stream(streamLength, totalDuration, dt, dx, boundary, DefaultParameters());
        SimulationResult result = stream.SimulateAdvectionAndDispersion();

        // The upstream boundary at t >= 10 s should be 0
        Assert.Equal(0.0, result.C[0, 50]);
    }
}
