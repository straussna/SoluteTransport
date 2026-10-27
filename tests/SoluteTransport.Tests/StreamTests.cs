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
    public void Constructor_ThrowsArgumentException_WhenASIsZeroWithExchange()
    {
        var p = DefaultParameters() with { AS = 0 };
        Assert.Throws<ArgumentException>(() =>
            new Stream(10, 10, 1, 1, [], p));
    }

    [Fact]
    public void Constructor_AcceptsZeroAS_WithoutExchange()
    {
        var p = DefaultParameters() with { AS = 0, Alpha = 0 };
        var result = new Stream(10, 10, 1, 1, [1.0], p).SimulateAdvectionAndDispersion();
        Assert.All(result.C.Cast<double>(), v => Assert.True(double.IsFinite(v)));
    }

    [Theory]
    [InlineData(1.0)]
    [InlineData(1.4)]
    public void Constructor_ThrowsArgumentException_WhenFewerThanTwoSegments(double streamLength)
    {
        Assert.Throws<ArgumentException>(() =>
            new Stream(streamLength, 10, 1, 1, [], DefaultParameters()));
    }

    [Fact]
    public void Constructor_ThrowsArgumentException_WhenDownstreamFluxWithoutDispersion()
    {
        var p = DefaultParameters() with { D = 0 };
        Assert.Throws<ArgumentException>(() =>
            new Stream(10, 10, 1, 1, [], p, downstreamDispersiveFlux: 0.1));
    }

    [Theory]
    [InlineData(-0.1, 0.0, 0.0)]
    [InlineData(0.5, -0.01, 0.0)]
    [InlineData(0.5, 0.0, -0.01)]
    public void Constructor_ThrowsArgumentException_WhenFlowIsNegative(double q, double qLin, double qLout)
    {
        var p = DefaultParameters() with { Q = q, QLin = qLin, QLout = qLout };
        Assert.Throws<ArgumentException>(() =>
            new Stream(10, 10, 1, 1, [], p));
    }

    [Fact]
    public void Constructor_ThrowsArgumentException_WhenOutflowDrainsReach()
    {
        // Q at the last segment centre is 0.5 - 0.06 * 9.5 < 0.
        var p = DefaultParameters() with { QLout = 0.06 };
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
        int numTimeLevels = (int)(totalDuration / dt) + 1;

        double[] boundary = new double[numTimeLevels];
        var stream = new Stream(streamLength, totalDuration, dt, dx, boundary, DefaultParameters());

        SimulationResult result = stream.SimulateAdvectionAndDispersion();

        Assert.Equal(numSegments, result.C.GetLength(0));
        Assert.Equal(numTimeLevels, result.C.GetLength(1));
        Assert.Equal(numSegments, result.CS.GetLength(0));
        Assert.Equal(numTimeLevels, result.CS.GetLength(1));
        Assert.Equal(numSegments, result.CSed.GetLength(0));
        Assert.Equal(numTimeLevels, result.CSed.GetLength(1));
    }

    [Fact]
    public void Simulate_AllZeroConcentrations_WhenNoInputAndNoTransport()
    {
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
        int lastTime = numTimeSteps;
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

        // Cell Peclet number Q dx / (A D) = 1, below the limit of 2 for
        // non-oscillating centred differences.
        var p = DefaultParameters() with { D = 0.5 };
        var stream = new Stream(streamLength, totalDuration, dt, dx, boundary, p);
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

        // Only supply 10 time levels of boundary data (less than numTimeSteps + 1)
        double[] shortBoundary = Enumerable.Repeat(1.0, 10).ToArray();
        double[] paddedBoundary = new double[numTimeSteps + 1];
        Array.Fill(paddedBoundary, 1.0, 0, 10);

        SimulationResult shortResult = new Stream(streamLength, totalDuration, dt, dx, shortBoundary, DefaultParameters())
            .SimulateAdvectionAndDispersion();
        SimulationResult paddedResult = new Stream(streamLength, totalDuration, dt, dx, paddedBoundary, DefaultParameters())
            .SimulateAdvectionAndDispersion();

        Assert.Equal(paddedResult.C.Cast<double>(), shortResult.C.Cast<double>());
    }

    [Fact]
    public void Simulate_SelectedTimeLevels_MatchFullRun()
    {
        var p = DefaultParameters() with { Kd = 0.2, LambdaHat = 0.01, Rho = 1.5 };
        var stream = new Stream(20, 200, 1, 0.5, [1.0, 1.0, 1.0, 0.5], p);
        int[] levels = [0, 7, 150, 200];

        SimulationResult full = stream.SimulateAdvectionAndDispersion();
        SimulationResult selected = stream.SimulateAdvectionAndDispersion(levels);

        Assert.Equal(levels, selected.TimeLevels);
        Assert.Equal(Enumerable.Range(0, 201), full.TimeLevels);
        Assert.Equal(levels.Length, selected.C.GetLength(1));
        for (int k = 0; k < levels.Length; k++)
        {
            for (int i = 0; i < stream.SegmentCount; i++)
            {
                Assert.Equal(full.C[i, levels[k]], selected.C[i, k]);
                Assert.Equal(full.CS[i, levels[k]], selected.CS[i, k]);
                Assert.Equal(full.CSed[i, levels[k]], selected.CSed[i, k]);
            }
        }
    }

    [Fact]
    public void Simulate_SelectedSegments_MatchFullRun()
    {
        var p = DefaultParameters() with { Kd = 0.2, LambdaHat = 0.01, Rho = 1.5 };
        var stream = new Stream(20, 200, 1, 0.5, [1.0, 1.0, 1.0, 0.5], p);
        int[] levels = [0, 7, 150, 200];
        int[] segments = [39, 0, 12, 12];

        SimulationResult full = stream.SimulateAdvectionAndDispersion(levels);
        SimulationResult selected = stream.SimulateAdvectionAndDispersion(levels, segments);

        Assert.Equal(segments, selected.Segments);
        Assert.Equal(Enumerable.Range(0, stream.SegmentCount), full.Segments);
        Assert.Equal(segments.Length, selected.C.GetLength(0));
        for (int s = 0; s < segments.Length; s++)
        {
            for (int k = 0; k < levels.Length; k++)
            {
                Assert.Equal(full.C[segments[s], k], selected.C[s, k]);
                Assert.Equal(full.CS[segments[s], k], selected.CS[s, k]);
                Assert.Equal(full.CSed[segments[s], k], selected.CSed[s, k]);
            }
        }
    }

    [Theory]
    [InlineData(-1)]
    [InlineData(40)]
    public void Simulate_ThrowsArgumentException_WhenSegmentOutOfRange(int segment)
    {
        var stream = new Stream(20, 10, 1, 0.5, [1.0], DefaultParameters());
        Assert.Throws<ArgumentException>(() => stream.SimulateAdvectionAndDispersion([10], [segment]));
    }

    [Fact]
    public void Simulate_RepeatedRuns_GiveIdenticalResults()
    {
        var stream = new Stream(10, 50, 1, 1, [1.0], DefaultParameters());
        Assert.Equal(
            stream.SimulateAdvectionAndDispersion().C.Cast<double>(),
            stream.SimulateAdvectionAndDispersion().C.Cast<double>());
    }

    [Theory]
    [InlineData(new[] { -1 })]
    [InlineData(new[] { 51 })]
    [InlineData(new[] { 5, 5 })]
    [InlineData(new[] { 10, 3 })]
    public void Simulate_ThrowsArgumentException_WhenTimeLevelsInvalid(int[] levels)
    {
        var stream = new Stream(10, 50, 1, 1, [], DefaultParameters());
        Assert.Throws<ArgumentException>(() => stream.SimulateAdvectionAndDispersion(levels));
    }
}
