namespace SoluteTransport;

/// <summary>
/// One-dimensional stream transport model that calculates solute advection and
/// dispersion using the Crank-Nicolson finite-difference scheme and the
/// Advection-Dispersion Equation (ADE) with transient storage and sediment
/// sorption.
/// </summary>
/// <remarks>
/// The numerical scheme uses the Thomas (tridiagonal matrix) algorithm to solve
/// the implicit finite-difference equations at each time step.  The upstream
/// boundary condition is a prescribed concentration time series; the downstream
/// boundary condition is zero concentration gradient (Neumann).
/// </remarks>
public class Stream
{
    // Grid and time-stepping
    private readonly ulong n;
    private readonly double dt;
    private readonly double dx;

    // Physical / chemical parameters
    private readonly double A;
    private readonly double AS;
    private readonly double D;
    private readonly double Q;
    private readonly double qLin;
    private readonly double alpha;
    private readonly double ChatS;
    private readonly double Kd;
    private readonly double lambda;
    private readonly double lambdaS;
    private readonly double lambdahat;
    private readonly double lambdahatS;
    private readonly double rho;
    private readonly double CL;

    // Concentration arrays  [segment, time-step]
    private readonly double[,] C;
    private readonly double[,] CS;
    private readonly double[,] Csed;

    // Pre-computed Crank-Nicolson tridiagonal coefficients:
    //   superDiagCoeff  – super-diagonal coefficient (downstream neighbour)
    //   subDiagCoeff    – sub-diagonal coefficient (upstream neighbour)
    //   mainDiagCoeff   – main diagonal coefficient (current segment)
    private readonly double gamma;
    private readonly double superDiagCoeff;
    private readonly double subDiagCoeff;
    private readonly double mainDiagCoeff;

    /// <summary>
    /// Initialises a new stream transport model.
    /// </summary>
    /// <param name="streamLength">Total length of the stream reach (m).</param>
    /// <param name="totalDuration">Total simulation duration (s).</param>
    /// <param name="dt">Time-step size (s).</param>
    /// <param name="dx">Spatial segment length (m).</param>
    /// <param name="upstreamBoundaryData">
    /// Upstream-boundary solute concentration (kg/m³) at each time step.
    /// Values beyond the array length are treated as zero.
    /// </param>
    /// <param name="parameters">Stream physical and chemical parameters.</param>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="upstreamBoundaryData"/> or
    /// <paramref name="parameters"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when <paramref name="dx"/> or <paramref name="dt"/> is
    /// not positive, or when <paramref name="parameters"/>.A or .AS is zero.
    /// </exception>
    public Stream(
        double streamLength,
        double totalDuration,
        double dt,
        double dx,
        double[] upstreamBoundaryData,
        StreamParameters parameters)
    {
        ArgumentNullException.ThrowIfNull(upstreamBoundaryData);
        ArgumentNullException.ThrowIfNull(parameters);
        if (dx <= 0) throw new ArgumentException("dx must be positive.", nameof(dx));
        if (dt <= 0) throw new ArgumentException("dt must be positive.", nameof(dt));
        if (parameters.A == 0) throw new ArgumentException("A must be non-zero.", nameof(parameters));
        if (parameters.AS == 0) throw new ArgumentException("AS must be non-zero.", nameof(parameters));

        ulong numberOfSegments = (ulong)(streamLength / dx);
        ulong numberOfTimeSteps = (ulong)(totalDuration / dt);

        this.n = numberOfSegments;
        this.dt = dt;
        this.dx = dx;

        this.A = parameters.A;
        this.AS = parameters.AS;
        this.CL = parameters.CL;
        this.D = parameters.D;
        this.Q = parameters.Q;
        this.qLin = parameters.QLin;
        this.alpha = parameters.Alpha;
        this.ChatS = parameters.ChatS;
        this.Kd = parameters.Kd;
        this.lambda = parameters.Lambda;
        this.lambdaS = parameters.LambdaS;
        this.lambdahat = parameters.LambdaHat;
        this.lambdahatS = parameters.LambdaHatS;
        this.rho = parameters.Rho;

        // Pre-compute Crank-Nicolson coefficients
        this.gamma = alpha * dt * A / AS;
        this.superDiagCoeff = dt / (2 * A * dx) * (Q / 2 - A * D / dx);
        this.subDiagCoeff = -dt / (2 * A * dx) * (Q / 2 + A * D / dx);
        this.mainDiagCoeff = 1 + dt / 2 * (
            (A * D + A * D) / (A * dx * dx)
            + qLin / A
            + alpha * (1 - alpha * dt * A / AS / (2 + alpha * dt * A / AS + dt * lambdahatS + dt * lambdaS))
            + rho * lambdahat * Kd * (1 - dt * lambdahat / (2 + dt * lambdahat))
            + lambda);

        this.C = new double[numberOfSegments, numberOfTimeSteps];
        this.CS = new double[numberOfSegments, numberOfTimeSteps];
        this.Csed = new double[numberOfSegments, numberOfTimeSteps];

        // Initialise interior segments to zero at t = 0
        for (int i = 1; i < C.GetLength(0); i++)
        {
            C[i, 0] = 0;
            CS[i, 0] = 0;
            Csed[i, 0] = 0;
        }

        // Apply upstream boundary condition
        for (int j = 0; j < C.GetLength(1); j++)
        {
            C[0, j] = j < upstreamBoundaryData.Length ? upstreamBoundaryData[j] : 0;
            CS[0, j] = 0;
            Csed[0, j] = 0;
        }
    }

    /// <summary>
    /// Runs the simulation and returns solute concentrations at every
    /// spatial segment and time step.
    /// </summary>
    /// <returns>
    /// A <see cref="SimulationResult"/> containing three matrices:
    /// main-channel concentration (C), storage-zone concentration (CS),
    /// and sediment-sorbed concentration (CSed).
    /// </returns>
    public SimulationResult SimulateAdvectionAndDispersion()
    {
        PropagateTime();
        return new SimulationResult(C, CS, Csed);
    }

    // ------------------------------------------------------------------
    // Private helpers
    // ------------------------------------------------------------------

    /// <summary>
    /// Right-hand-side contribution from segment <paramref name="i"/> at
    /// time step <paramref name="j"/>, incorporating storage-zone and
    /// sediment exchange.
    /// </summary>
    private double R(int i, int j)
    {
        return C[i, j]
            + dt / 2 * (
                superDiagCoeff + qLin / A * CL
                + alpha * ((2 - gamma - dt * lambdahatS - dt * lambdaS) * CS[i, j]
                           + gamma * C[i, j]
                           + 2 * dt * lambdahatS * ChatS)
                  / (2 + gamma + dt * lambdahatS + dt * lambdaS))
            + rho * lambdahat
              * ((2 - dt * lambdahat) * Csed[i, j]
                 + dt * lambdahat * Kd * C[i, j])
              / (2 + dt * lambdahat);
    }

    /// <summary>Advances the simulation over all time steps.</summary>
    private void PropagateTime()
    {
        for (int j = 0; j < C.GetLength(1) - 1; j++)
        {
            IncrementTimeForward(j);
        }
    }

    /// <summary>
    /// Propagates C, CS, and Csed from time step <paramref name="j"/> to
    /// <paramref name="j"/> + 1 using the Thomas (tridiagonal matrix)
    /// algorithm.
    /// See https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    /// </summary>
    private void IncrementTimeForward(int j)
    {
        int m = C.GetLength(0) - 1; // interior segment count

        // RHS vector
        double[] d = new double[m];
        // Tridiagonal matrix bands: sub-diagonal (a), main diagonal (b), super-diagonal (c)
        double[] a = new double[m];
        double[] b = new double[m];
        double[] c = new double[m];

        for (int i = 0; i < m; i++)
        {
            d[i] = R(i, j);
            a[i] = subDiagCoeff;
            b[i] = mainDiagCoeff;
            c[i] = superDiagCoeff;
        }

        // Upstream boundary condition: prescribed concentration
        d[0] = R(1, j) - subDiagCoeff * C[0, j + 1];

        // Downstream boundary condition: zero concentration gradient
        b[m - 1] = mainDiagCoeff + superDiagCoeff;

        // Thomas algorithm – forward sweep
        double[] cp = new double[m - 1];
        double[] dp = new double[m];

        cp[0] = c[0] / b[0];
        dp[0] = d[0] / b[0];

        for (int i = 1; i < m - 1; i++)
        {
            double denom = b[i] - a[i] * cp[i - 1];
            cp[i] = c[i] / denom;
            dp[i] = (d[i] - a[i] * dp[i - 1]) / denom;
        }
        dp[m - 1] = (d[m - 1] - a[m - 1] * dp[m - 2]) / (b[m - 1] - a[m - 1] * cp[m - 2]);

        // Back-substitution
        C[C.GetLength(0) - 1, j + 1] = dp[m - 1];
        UpdateStorageAndSediment(C.GetLength(0) - 1, j);

        for (int i = C.GetLength(0) - 2; i >= 1; i--)
        {
            C[i, j + 1] = dp[i - 1] - cp[i - 1] * C[i + 1, j + 1];
            UpdateStorageAndSediment(i, j);
        }
    }

    /// <summary>
    /// Updates the transient-storage and sediment concentrations at segment
    /// <paramref name="i"/> from time step <paramref name="j"/> to j + 1.
    /// </summary>
    private void UpdateStorageAndSediment(int i, int j)
    {
        CS[i, j + 1] = ((2 - gamma - dt * lambdahatS - dt * lambdaS) * CS[i, j]
                        + gamma * C[i, j]
                        + gamma * C[i, j + 1]
                        + 2 * dt * lambdahatS * ChatS)
                       / (2 + gamma + dt * lambdahatS + dt * lambdaS);

        Csed[i, j + 1] = ((2 - dt * lambdahat) * Csed[i, j]
                          + dt * lambdahat * Kd * (C[i, j] + C[i, j + 1]))
                         / (2 + dt * lambdahat);
    }
}
