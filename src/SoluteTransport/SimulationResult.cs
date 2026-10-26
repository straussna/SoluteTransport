namespace SoluteTransport;

/// <summary>
/// Concentrations at every segment and time level produced by
/// <see cref="Stream.SimulateAdvectionAndDispersion"/>.
/// </summary>
/// <param name="C">
/// Main-channel solute concentration (kg/m³).
/// First index: segment <c>i</c>, centred at <c>(i + 0.5) * dx</c> from the
/// upstream boundary; second index: time level <c>j</c>, at time
/// <c>j * dt</c>, from 0 to the total duration inclusive.
/// </param>
/// <param name="CS">
/// Transient-storage zone solute concentration (kg/m³).
/// Indices identical to <paramref name="C"/>.
/// </param>
/// <param name="CSed">
/// Sediment-sorbed solute concentration (kg/kg).
/// Indices identical to <paramref name="C"/>.
/// </param>
public record SimulationResult(double[,] C, double[,] CS, double[,] CSed);
