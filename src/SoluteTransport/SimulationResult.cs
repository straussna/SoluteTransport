namespace SoluteTransport;

/// <summary>
/// Concentrations at every spatial segment and time step produced by
/// <see cref="Stream.SimulateAdvectionAndDispersion"/>.
/// </summary>
/// <param name="C">
/// Main-channel solute concentration (kg/m³).
/// First index: spatial segment (0 = upstream); second index: time step.
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
