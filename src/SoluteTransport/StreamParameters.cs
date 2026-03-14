namespace SoluteTransport;

/// <summary>
/// Physical and chemical parameters that describe a stream reach for the 1-D
/// Advection-Dispersion Equation (ADE) with transient storage and sediment
/// sorption.
/// </summary>
/// <param name="A">Main-channel cross-sectional area (m²).</param>
/// <param name="AS">Transient-storage zone cross-sectional area (m²).</param>
/// <param name="CL">Lateral inflow solute concentration (kg/m³).</param>
/// <param name="D">Longitudinal dispersion coefficient (m²/s).</param>
/// <param name="Q">Volumetric flow rate (m³/s).</param>
/// <param name="QLin">Lateral inflow rate per unit length (m²/s).</param>
/// <param name="Alpha">Exchange rate between main channel and storage zone (1/s).</param>
/// <param name="ChatS">Equilibrium solute concentration in the storage zone (kg/m³).</param>
/// <param name="Kd">Sediment-water partition coefficient (m³/kg).</param>
/// <param name="Lambda">First-order decay rate in main channel (1/s).</param>
/// <param name="LambdaS">First-order decay rate in storage zone (1/s).</param>
/// <param name="LambdaHat">First-order desorption rate from sediment (1/s).</param>
/// <param name="LambdaHatS">First-order desorption rate in storage zone (1/s).</param>
/// <param name="Rho">Sediment bulk density (kg/m³).</param>
public record StreamParameters(
    double A,
    double AS,
    double CL,
    double D,
    double Q,
    double QLin,
    double Alpha,
    double ChatS,
    double Kd,
    double Lambda,
    double LambdaS,
    double LambdaHat,
    double LambdaHatS,
    double Rho);
