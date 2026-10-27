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
/// The tridiagonal coefficients are constant in time, so the elimination
/// factors are computed once; the simulation holds two time levels in memory
/// and stores only the requested output levels.
/// </para>
/// <para>
/// Flow is steady and varies along the reach by the water balance
/// <c>dQ/dx = qLIN - qLOUT</c> (Runkel and Broshears 1991, equation A20).
/// The flow rate at each segment centre follows Runkel and Broshears (1991),
/// equation (31), starting from <c>Q</c> at the upstream boundary. Lateral
/// outflow removes water at the main-channel concentration, so it changes
/// the solute balance only through the flow rate.
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
    private readonly double[] upstreamBoundary;

    // Right-hand side of equation (17) at segment i:
    //   explicitLower[i] C[i-1] + explicitCentre[i] C[i] + explicitUpper[i] C[i+1]
    //   + storageWeight CS[i] + sedimentWeight Csed[i] + sourceTerm,
    // plus upstreamTerm (Cbc[j] + Cbc[j+1]) in the first segment and
    // downstreamTerm in the last.
    private readonly double[] explicitLower;
    private readonly double[] explicitCentre;
    private readonly double[] explicitUpper;
    private readonly double storageWeight;
    private readonly double sedimentWeight;
    private readonly double sourceTerm;
    private readonly double upstreamTerm;
    private readonly double downstreamTerm;

    // Equations (25) and (26) as
    //   CS[j+1]   = storageOld CS[j]     + storageFromC (C[j] + C[j+1]) + storageConstant
    //   Csed[j+1] = sedimentOld Csed[j]  + sedimentFromC (C[j] + C[j+1]).
    private readonly double storageOld;
    private readonly double storageFromC;
    private readonly double storageConstant;
    private readonly double sedimentOld;
    private readonly double sedimentFromC;

    // Thomas-algorithm factors of the constant Crank-Nicolson matrix:
    // sub-diagonal, modified super-diagonal and reciprocal pivots.
    private readonly double[] subDiagonal;
    private readonly double[] modifiedUpper;
    private readonly double[] reciprocalPivot;

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
    /// less than one time step, when A is not positive, when D, QLin or QLout
    /// is negative, when the flow rate is negative at the upstream boundary or
    /// any segment centre, when AS is not positive while Alpha is non-zero, or
    /// when D is zero while <paramref name="downstreamDispersiveFlux"/> is
    /// non-zero.
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
        if (parameters.Q < 0) throw new ArgumentException("Q must not be negative.", nameof(parameters));
        if (parameters.QLin < 0) throw new ArgumentException("QLin must not be negative.", nameof(parameters));
        if (parameters.QLout < 0) throw new ArgumentException("QLout must not be negative.", nameof(parameters));
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

        StreamParameters p = parameters;
        upstreamBoundary = upstreamBoundaryData;

        // Flow rate at segment centres, Runkel and Broshears (1991) equation (31),
        // with equal segment lengths and uniform lateral flows.
        double lateralNet = p.QLin - p.QLout;
        double[] flow = new double[segmentCount];
        for (int i = 0; i < segmentCount; i++)
        {
            flow[i] = p.Q + lateralNet * (i + 0.5) * dx;
        }
        if (flow[segmentCount - 1] < 0)
            throw new ArgumentException("QLout must not drain the flow rate below zero within the reach.", nameof(parameters));

        // Interior segments: Runkel (1998) equation (12).
        double[] adv = new double[segmentCount];
        double disp = p.D / (dx * dx);
        // Transport operator L[C] = lower*C[i-1] + diag*C[i] + upper*C[i+1].
        double[] lower = new double[segmentCount];
        double[] diag = new double[segmentCount];
        double[] upper = new double[segmentCount];
        for (int i = 0; i < segmentCount; i++)
        {
            adv[i] = flow[i] / (2 * p.A * dx);
            lower[i] = adv[i] + disp;
            diag[i] = -2 * disp;
            upper[i] = -adv[i] + disp;
        }

        // First segment: equations (B11) and (B14).
        //   advection  -(Q[0] / A dx) ((C[0] + C[1]) / 2 - Cbc)
        //   dispersion (D / dx²) (C[1] - 3 C[0] + 2 Cbc)
        lower[0] = 0;
        diag[0] = -adv[0] - 3 * disp;
        upper[0] = -adv[0] + disp;
        double upstreamWeight = 2 * adv[0] + 2 * disp;

        // Last segment: equations (B16)–(B18), with C[N] = C[N-1] + dx DSBOUND / D.
        //   advection  -(Q[N-1] / 2 A dx) (C[N-1] - C[N-2] + dx DSBOUND / D)
        //   dispersion DSBOUND / dx - (D / dx²) (C[N-1] - C[N-2])
        int last = segmentCount - 1;
        diag[last] = -adv[last] - disp;
        upper[last] = 0;
        double downstreamSource = downstreamDispersiveFlux == 0
            ? 0
            : -flow[last] * downstreamDispersiveFlux / (2 * p.A * p.D) + downstreamDispersiveFlux / dx;

        // Decoupling constants, Runkel (1998) equations (25)–(27).
        double gamma = p.Alpha == 0 ? 0 : p.Alpha * dt * p.A / p.AS;
        double storageDenominator = 2 + gamma + dt * p.LambdaHatS + dt * p.LambdaS;
        double sedimentDenominator = 2 + dt * p.LambdaHat;
        storageOld = (2 - gamma - dt * p.LambdaHatS - dt * p.LambdaS) / storageDenominator;
        storageFromC = gamma / storageDenominator;
        storageConstant = 2 * dt * p.LambdaHatS * p.ChatS / storageDenominator;
        sedimentOld = (2 - dt * p.LambdaHat) / sedimentDenominator;
        sedimentFromC = dt * p.LambdaHat * p.Kd / sedimentDenominator;

        // Reaction part of F' (equation 29).
        double lateral = p.QLin / p.A;
        double sorption = p.Rho * p.LambdaHat;
        double reaction = lateral
            + p.Alpha * (1 - storageFromC)
            + sorption * (p.Kd - sedimentFromC)
            + p.Lambda;

        // Equation (17): the implicit half on the left, the explicit half
        // G[C, CS, Csed] (equation 15) and the known parts of equation (30)
        // on the right.
        double half = dt / 2;
        storageWeight = half * p.Alpha * (1 + storageOld);
        sedimentWeight = half * sorption * (1 + sedimentOld);
        sourceTerm = half * (2 * lateral * p.CL + p.Alpha * storageConstant);
        upstreamTerm = half * upstreamWeight;
        downstreamTerm = dt * downstreamSource;

        explicitLower = new double[segmentCount];
        explicitCentre = new double[segmentCount];
        explicitUpper = new double[segmentCount];
        subDiagonal = new double[segmentCount];
        modifiedUpper = new double[segmentCount];
        reciprocalPivot = new double[segmentCount];
        for (int i = 0; i < segmentCount; i++)
        {
            explicitLower[i] = half * lower[i];
            explicitCentre[i] = 1 - half * (reaction - diag[i]);
            explicitUpper[i] = half * upper[i];

            double a = -half * lower[i];
            double b = 1 + half * (reaction - diag[i]);
            double c = -half * upper[i];
            double pivot = i == 0 ? b : b - a * modifiedUpper[i - 1];
            subDiagonal[i] = a;
            reciprocalPivot[i] = 1 / pivot;
            modifiedUpper[i] = c / pivot;
        }
    }

    /// <summary>Number of segments in the reach.</summary>
    public int SegmentCount => segmentCount;

    /// <summary>
    /// Number of time steps; time levels run from 0 to
    /// <see cref="TimeStepCount"/> inclusive.
    /// </summary>
    public int TimeStepCount => timeStepCount;

    /// <summary>
    /// Runs the simulation and returns solute concentrations at every
    /// segment and time level.
    /// </summary>
    /// <returns>
    /// A <see cref="SimulationResult"/> containing main-channel (C),
    /// storage-zone (CS) and streambed-sediment (CSed) concentrations at time
    /// levels 0 to <see cref="TimeStepCount"/>.
    /// </returns>
    public SimulationResult SimulateAdvectionAndDispersion() =>
        SimulateAdvectionAndDispersion(Enumerable.Range(0, timeStepCount + 1).ToArray());

    /// <summary>
    /// Runs the simulation up to the last requested time level and returns
    /// solute concentrations at every segment for the requested levels only.
    /// </summary>
    /// <param name="timeLevels">
    /// Strictly increasing time levels <c>j</c> (time <c>j * dt</c>), each
    /// between 0 and <see cref="TimeStepCount"/> inclusive.
    /// </param>
    /// <returns>
    /// A <see cref="SimulationResult"/> whose second index <c>k</c> is the
    /// time level <c>timeLevels[k]</c>.
    /// </returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="timeLevels"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when a time level is out of range or the levels are not
    /// strictly increasing.
    /// </exception>
    public SimulationResult SimulateAdvectionAndDispersion(IReadOnlyList<int> timeLevels) =>
        SimulateAdvectionAndDispersion(timeLevels, Enumerable.Range(0, segmentCount).ToArray());

    /// <summary>
    /// Runs the simulation up to the last requested time level and returns
    /// solute concentrations at the requested segments and time levels only.
    /// </summary>
    /// <param name="timeLevels">
    /// Strictly increasing time levels <c>j</c> (time <c>j * dt</c>), each
    /// between 0 and <see cref="TimeStepCount"/> inclusive.
    /// </param>
    /// <param name="segments">
    /// Segment indices, each between 0 and <see cref="SegmentCount"/> - 1, in
    /// any order.
    /// </param>
    /// <returns>
    /// A <see cref="SimulationResult"/> whose first index <c>s</c> is the
    /// segment <c>segments[s]</c> and second index <c>k</c> the time level
    /// <c>timeLevels[k]</c>.
    /// </returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="timeLevels"/> or <paramref name="segments"/>
    /// is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when a time level or segment is out of range or the levels are
    /// not strictly increasing.
    /// </exception>
    public SimulationResult SimulateAdvectionAndDispersion(IReadOnlyList<int> timeLevels, IReadOnlyList<int> segments)
    {
        ArgumentNullException.ThrowIfNull(timeLevels);
        ArgumentNullException.ThrowIfNull(segments);
        for (int k = 0; k < timeLevels.Count; k++)
        {
            if (timeLevels[k] < 0 || timeLevels[k] > timeStepCount)
                throw new ArgumentException($"Time level {timeLevels[k]} is outside 0 to {timeStepCount}.", nameof(timeLevels));
            if (k > 0 && timeLevels[k] <= timeLevels[k - 1])
                throw new ArgumentException("Time levels must be strictly increasing.", nameof(timeLevels));
        }
        foreach (int segment in segments)
        {
            if (segment < 0 || segment >= segmentCount)
                throw new ArgumentException($"Segment {segment} is outside 0 to {segmentCount - 1}.", nameof(segments));
        }

        int[] levels = timeLevels.ToArray();
        int[] outputSegments = segments.ToArray();
        var resultC = new double[outputSegments.Length, levels.Length];
        var resultCS = new double[outputSegments.Length, levels.Length];
        var resultCSed = new double[outputSegments.Length, levels.Length];

        var current = new Level(segmentCount);
        var next = new Level(segmentCount);

        int output = 0;
        for (int j = 0; output < levels.Length; j++)
        {
            if (j > 0)
            {
                IncrementTimeForward(j - 1, current, next);
                (current, next) = (next, current);
            }
            if (levels[output] == j)
            {
                for (int s = 0; s < outputSegments.Length; s++)
                {
                    int i = outputSegments[s];
                    resultC[s, output] = current.C[i];
                    resultCS[s, output] = current.CS[i];
                    resultCSed[s, output] = current.Csed[i];
                }
                output++;
            }
        }

        return new SimulationResult(resultC, resultCS, resultCSed, outputSegments, levels);
    }

    /// <summary>Main-channel, storage-zone and sediment concentrations at one time level.</summary>
    private sealed class Level(int segmentCount)
    {
        public readonly double[] C = new double[segmentCount];
        public readonly double[] CS = new double[segmentCount];
        public readonly double[] Csed = new double[segmentCount];
    }

    private double UpstreamConcentration(int j) =>
        j < upstreamBoundary.Length ? upstreamBoundary[j] : 0;

    /// <summary>
    /// Advances <paramref name="current"/> at time level <paramref name="j"/>
    /// to <paramref name="next"/> at <paramref name="j"/> + 1. The forward
    /// sweep of the Thomas algorithm builds the right-hand side of equation
    /// (17) as it goes; the backward sweep finishes each C[i] and applies
    /// equations (25) and (26) to the same segment.
    /// </summary>
    private void IncrementTimeForward(int j, Level current, Level next)
    {
        int n = segmentCount;
        int last = n - 1;
        double[] c = current.C, cs = current.CS, sed = current.Csed;
        double[] x = next.C, nextCS = next.CS, nextSed = next.Csed;
        double[] el = explicitLower, ec = explicitCentre, eu = explicitUpper;
        double[] sub = subDiagonal, mu = modifiedUpper, rp = reciprocalPivot;
        double sw = storageWeight, dw = sedimentWeight, source = sourceTerm;

        double rhs = ec[0] * c[0] + eu[0] * c[1] + sw * cs[0] + dw * sed[0] + source
                     + upstreamTerm * (UpstreamConcentration(j) + UpstreamConcentration(j + 1));
        x[0] = rhs * rp[0];
        for (int i = 1; i < last; i++)
        {
            rhs = el[i] * c[i - 1] + ec[i] * c[i] + eu[i] * c[i + 1]
                  + sw * cs[i] + dw * sed[i] + source;
            x[i] = (rhs - sub[i] * x[i - 1]) * rp[i];
        }
        rhs = el[last] * c[last - 1] + ec[last] * c[last]
              + sw * cs[last] + dw * sed[last] + source + downstreamTerm;
        x[last] = (rhs - sub[last] * x[last - 1]) * rp[last];

        double so = storageOld, sc = storageFromC, sk = storageConstant;
        double ho = sedimentOld, hc = sedimentFromC;
        for (int i = last; i >= 0; i--)
        {
            if (i < last) x[i] -= mu[i] * x[i + 1];
            double sum = c[i] + x[i];
            nextCS[i] = so * cs[i] + sc * sum + sk;
            nextSed[i] = ho * sed[i] + hc * sum;
        }
    }
}
