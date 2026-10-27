namespace SoluteTransport;

/// <summary>
/// Concentrations at the output segments and time levels of a
/// <see cref="Stream"/> simulation.
/// </summary>
/// <param name="C">
/// Main-channel solute concentration (kg/m³).
/// First index: output <c>s</c>, at segment <c>Segments[s]</c>, centred at
/// <c>(Segments[s] + 0.5) * dx</c> from the upstream boundary; second index:
/// output <c>k</c>, at time level <c>TimeLevels[k]</c> (time
/// <c>TimeLevels[k] * dt</c>).
/// </param>
/// <param name="CS">
/// Transient-storage zone solute concentration (kg/m³).
/// Indices identical to <paramref name="C"/>.
/// </param>
/// <param name="CSed">
/// Sediment-sorbed solute concentration (kg/kg).
/// Indices identical to <paramref name="C"/>.
/// </param>
/// <param name="Segments">Segment index of each output.</param>
/// <param name="TimeLevels">Time level of each output, strictly increasing.</param>
public record SimulationResult(double[,] C, double[,] CS, double[,] CSed, int[] Segments, int[] TimeLevels);
