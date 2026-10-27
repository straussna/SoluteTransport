namespace SoluteTransport;

/// <summary>
/// Physical and chemical parameters that describe a stream reach for the 1-D
/// Advection-Dispersion Equation (ADE) with transient storage and sediment
/// sorption. Symbols and meanings follow Runkel (1998), equations (1)–(5).
/// </summary>
/// <param name="A">Main-channel cross-sectional area (m²).</param>
/// <param name="AS">Storage-zone cross-sectional area (m²).</param>
/// <param name="CL">Lateral inflow solute concentration (kg/m³).</param>
/// <param name="D">Longitudinal dispersion coefficient (m²/s).</param>
/// <param name="Q">Volumetric flow rate at the upstream boundary, QSTART (m³/s).</param>
/// <param name="QLin">Lateral inflow rate per unit length (m²/s).</param>
/// <param name="Alpha">Storage-zone exchange coefficient α (1/s).</param>
/// <param name="ChatS">Background storage-zone solute concentration Ĉ_S (kg/m³).</param>
/// <param name="Kd">Distribution coefficient K_d (m³/kg).</param>
/// <param name="Lambda">Main-channel first-order decay coefficient λ (1/s).</param>
/// <param name="LambdaS">Storage-zone first-order decay coefficient λ_S (1/s).</param>
/// <param name="LambdaHat">Main-channel sorption rate coefficient λ̂ (1/s).</param>
/// <param name="LambdaHatS">Storage-zone sorption rate coefficient λ̂_S (1/s).</param>
/// <param name="Rho">Mass of accessible sediment per volume of water ρ (kg/m³).</param>
/// <param name="QLout">Lateral outflow rate per unit length (m²/s).</param>
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
    double Rho,
    double QLout = 0);
