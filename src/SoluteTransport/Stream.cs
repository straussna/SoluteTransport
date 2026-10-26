namespace SoluteTransport;

/// <summary>
/// One-dimensional stream transport model (OTIS, Runkel 1998) that solves the
/// Advection-Dispersion Equation with lateral inflow, transient storage,
/// streambed sorption and first-order decay by the Crank-Nicolson
/// finite-difference method.
/// </summary>
/// <remarks>
/// <para>
/// The reach is divided into segments of equal length <c>dx</c>; segment
/// <c>i</c> (0-based) is centred at <c>(i + 0.5) * dx</c>. Main-channel,
/// storage-zone and streambed-sediment equations are decoupled as in
/// Runkel (1998) equations (25)–(30), so each time step is one tridiagonal
/// solve (Thomas algorithm) followed by explicit storage and sediment updates.
/// </para>
/// <para>
/// Boundary conditions follow Runkel and Broshears (1991), appendix B: the
/// upstream boundary is a fixed concentration <c>Cbc</c> at the upstream face
/// of the first segment (equations B10–B14), and the downstream boundary is a
/// fixed dispersive flux <c>D ∂C/∂x</c> at the downstream face of the last
/// segment (equations B15–B18); a flux of zero is a zero-gradient condition.
/// </para>
/// </remarks>
public class Stream
{
    private readonly int segmentCount;
    private readonly int timeStepCount;
    private readonly double dt;
    private readonly StreamParameters p;
    private readonly double[] upstreamBoundary;

    // Concentration arrays [segment, time level]
    private readonly double[,] C;
    private readonly double[,] CS;
    private readonly double[,] Csed;

    // Transport operator L[C] = lower*C[i-1] + diag*C[i] + upper*C[i+1],
    // with segment-specific entries at the boundaries.
    private readonly double[] lower;
    private readonly double[] diag;
    private readonly double[] upper;

    // Weight of Cbc in L[C] for the first segment, and the constant
    // downstream-flux source in L[C] for the last segment.
    private readonly double upstreamWeight;
    private readonly double downstreamSource;

    // Decoupling constants, Runkel (1998) equations (25)–(27).
    private readonly double gamma;
    private readonly double storageDenominator;
    private readonly double sedimentDenominator;

    /// <summary>
    /// Initialises a new stream transport model.
    /// </summary>
    /// <param name="streamLength">Total length of the stream reach (m); rounded to a whole number of segments.</param>
    /// <param name="totalDuration">Total simulation duration (s); rounded to a whole number of time steps.</param>
    /// <param name="dt">Time-step size (s).</param>
    /// <param name="dx">Segment length (m).</param>
    /// <param name="upstreamBoundaryData">
    /// Upstream-boundary solute concentration <c>Cbc</c> (kg/m³) at each time
    /// level <c>j</c> (time <c>j * dt</c>). Values beyond the array length are
    /// treated as zero.
    /// </param>
    /// <param name="parameters">Stream physical and chemical parameters.</param>
    /// <param name="downstreamDispersiveFlux">
    /// Dispersive flux <c>D ∂C/∂x</c> (kg/m²/s) at the downstream boundary
    /// (DSBOUND); zero gives a zero-gradient boundary.
    /// </param>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="upstreamBoundaryData"/> or
    /// <paramref name="parameters"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when <paramref name="dx"/> or <paramref name="dt"/> is not
    /// positive, when the reach holds fewer than two segments or the duration
    /// less than one time step, when A is not positive, when D is negative,
    /// when AS is not positive while Alpha is non-zero, or when D is zero while
    /// <paramref name="downstreamDispersiveFlux"/> is non-zero.
    /// </exception>
    public Stream(
        double streamLength,
        double totalDuration,
        double dt,
        double dx,
        double[] upstreamBoundaryData,
        StreamParameters parameters,
        double downstreamDispersiveFlux = 0)
    {
        ArgumentNullException.ThrowIfNull(upstreamBoundaryData);
        ArgumentNullException.ThrowIfNull(parameters);
        if (dx <= 0) throw new ArgumentException("dx must be positive.", nameof(dx));
        if (dt <= 0) throw new ArgumentException("dt must be positive.", nameof(dt));
        if (parameters.A <= 0) throw new ArgumentException("A must be positive.", nameof(parameters));
        if (parameters.D < 0) throw new ArgumentException("D must not be negative.", nameof(parameters));
        if (parameters.Alpha != 0 && parameters.AS <= 0)
            throw new ArgumentException("AS must be positive when Alpha is non-zero.", nameof(parameters));
        if (downstreamDispersiveFlux != 0 && parameters.D == 0)
            throw new ArgumentException("D must be positive when the downstream dispersive flux is non-zero.", nameof(downstreamDispersiveFlux));

        segmentCount = (int)Math.Round(streamLength / dx);
        timeStepCount = (int)Math.Round(totalDuration / dt);
        if (segmentCount < 2)
            throw new ArgumentException("streamLength must hold at least two segments of length dx.", nameof(streamLength));
        if (timeStepCount < 1)
            throw new ArgumentException("totalDuration must be at least one time step.", nameof(totalDuration));

        this.dt = dt;
        this.p = parameters;
        this.upstreamBoundary = upstreamBoundaryData;

        C = new double[segmentCount, timeStepCount + 1];
        CS = new double[segmentCount, timeStepCount + 1];
        Csed = new double[segmentCount, timeStepCount + 1];

        // Interior segments: Runkel (1998) equation (12).
        double adv = p.Q / (2 * p.A * dx);
        double disp = p.D / (dx * dx);
        lower = new double[segmentCount];
        diag = new double[segmentCount];
        upper = new double[segmentCount];
        for (int i = 0; i < segmentCount; i++)
        {
            lower[i] = adv + disp;
            diag[i] = -2 * disp;
            upper[i] = -adv + disp;
        }

        // First segment: equations (B11) and (B14).
        //   advection  -(Q / A dx) ((C[0] + C[1]) / 2 - Cbc)
        //   dispersion (D / dx²) (C[1] - 3 C[0] + 2 Cbc)
        lower[0] = 0;
        diag[0] = -adv - 3 * disp;
        upper[0] = -adv + disp;
        upstreamWeight = 2 * adv + 2 * disp;

        // Last segment: equations (B16)–(B18), with C[N] = C[N-1] + dx DSBOUND / D.
        //   advection  -(Q / 2 A dx) (C[N-1] - C[N-2] + dx DSBOUND / D)
        //   dispersion DSBOUND / dx - (D / dx²) (C[N-1] - C[N-2])
        int last = segmentCount - 1;
        diag[last] = -adv - disp;
        upper[last] = 0;
        downstreamSource = downstreamDispersiveFlux == 0
            ? 0
            : -p.Q * downstreamDispersiveFlux / (2 * p.A * p.D) + downstreamDispersiveFlux / dx;

        gamma = p.Alpha == 0 ? 0 : p.Alpha * dt * p.A / p.AS;
        storageDenominator = 2 + gamma + dt * p.LambdaHatS + dt * p.LambdaS;
        sedimentDenominator = 2 + dt * p.LambdaHat;
    }

    /// <summary>
    /// Runs the simulation and returns solute concentrations at every
    /// segment and time level.
    /// </summary>
    /// <returns>
    /// A <see cref="SimulationResult"/> containing main-channel (C),
    /// storage-zone (CS) and streambed-sediment (CSed) concentrations.
    /// </returns>
    public SimulationResult SimulateAdvectionAndDispersion()
    {
        for (int j = 0; j < timeStepCount; j++)
        {
            IncrementTimeForward(j);
        }
        return new SimulationResult(C, CS, Csed);
    }

    private double UpstreamConcentration(int j) =>
        j < upstreamBoundary.Length ? upstreamBoundary[j] : 0;

    /// <summary>
    /// G[C, CS, Csed] at segment <paramref name="i"/> and time level
    /// <paramref name="j"/> (Runkel 1998, equation (15)), including the
    /// boundary terms of the transport operator.
    /// </summary>
    private double G(int i, int j)
    {
        double transport = diag[i] * C[i, j];
        if (i > 0) transport += lower[i] * C[i - 1, j];
        if (i < segmentCount - 1) transport += upper[i] * C[i + 1, j];
        if (i == 0) transport += upstreamWeight * UpstreamConcentration(j);
        if (i == segmentCount - 1) transport += downstreamSource;

        return transport
            + p.QLin / p.A * (p.CL - C[i, j])
            + p.Alpha * (CS[i, j] - C[i, j])
            + p.Rho * p.LambdaHat * (Csed[i, j] - p.Kd * C[i, j])
            - p.Lambda * C[i, j];
    }

    /// <summary>
    /// Advances C, CS and Csed from time level <paramref name="j"/> to
    /// <paramref name="j"/> + 1: solves equation (17) with the decoupled
    /// coefficients of equations (18), (20), (29) and (30), then applies
    /// equations (25) and (26).
    /// </summary>
    private void IncrementTimeForward(int j)
    {
        int n = segmentCount;
        double half = dt / 2;

        // Reaction part of F' (equation 29).
        double reaction = p.QLin / p.A
            + p.Alpha * (1 - gamma / storageDenominator)
            + p.Rho * p.LambdaHat * p.Kd * (1 - dt * p.LambdaHat / sedimentDenominator)
            + p.Lambda;

        double[] a = new double[n];
        double[] b = new double[n];
        double[] c = new double[n];
        double[] r = new double[n];

        for (int i = 0; i < n; i++)
        {
            a[i] = -half * lower[i];
            b[i] = 1 + half * (reaction - diag[i]);
            c[i] = -half * upper[i];

            // Known parts of CS and Csed at j + 1, equation (30).
            double storageKnown = ((2 - gamma - dt * p.LambdaHatS - dt * p.LambdaS) * CS[i, j]
                                   + gamma * C[i, j]
                                   + 2 * dt * p.LambdaHatS * p.ChatS)
                                  / storageDenominator;
            double sedimentKnown = ((2 - dt * p.LambdaHat) * Csed[i, j]
                                    + dt * p.LambdaHat * p.Kd * C[i, j])
                                   / sedimentDenominator;

            r[i] = C[i, j] + half * (
                G(i, j)
                + p.QLin / p.A * p.CL
                + p.Alpha * storageKnown
                + p.Rho * p.LambdaHat * sedimentKnown);
        }

        // Boundary terms at the advanced time level.
        r[0] += half * upstreamWeight * UpstreamConcentration(j + 1);
        r[n - 1] += half * downstreamSource;

        double[] next = SolveTridiagonal(a, b, c, r);

        for (int i = 0; i < n; i++)
        {
            C[i, j + 1] = next[i];

            CS[i, j + 1] = ((2 - gamma - dt * p.LambdaHatS - dt * p.LambdaS) * CS[i, j]
                            + gamma * C[i, j]
                            + gamma * C[i, j + 1]
                            + 2 * dt * p.LambdaHatS * p.ChatS)
                           / storageDenominator;

            Csed[i, j + 1] = ((2 - dt * p.LambdaHat) * Csed[i, j]
                              + dt * p.LambdaHat * p.Kd * (C[i, j] + C[i, j + 1]))
                             / sedimentDenominator;
        }
    }

    /// <summary>
    /// Solves a tridiagonal system by the Thomas algorithm.
    /// See https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    /// </summary>
    /// <param name="a">Sub-diagonal; <c>a[0]</c> is unused.</param>
    /// <param name="b">Main diagonal.</param>
    /// <param name="c">Super-diagonal; the last entry is unused.</param>
    /// <param name="d">Right-hand side.</param>
    private static double[] SolveTridiagonal(double[] a, double[] b, double[] c, double[] d)
    {
        int n = d.Length;
        double[] cp = new double[n];
        double[] dp = new double[n];

        cp[0] = c[0] / b[0];
        dp[0] = d[0] / b[0];
        for (int i = 1; i < n; i++)
        {
            double denom = b[i] - a[i] * cp[i - 1];
            cp[i] = c[i] / denom;
            dp[i] = (d[i] - a[i] * dp[i - 1]) / denom;
        }

        double[] x = new double[n];
        x[n - 1] = dp[n - 1];
        for (int i = n - 2; i >= 0; i--)
        {
            x[i] = dp[i] - cp[i] * x[i + 1];
        }
        return x;
    }
}
