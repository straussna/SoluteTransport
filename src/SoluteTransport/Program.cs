using SoluteTransport;

// ---------------------------------------------------------------------------
// Example: simulate a conservative tracer pulse injected at the upstream end
// of a 50-metre stream reach for 1 hour.
// ---------------------------------------------------------------------------

const double streamLength = 50.0;   // m
const double totalDuration = 3600.0; // s
const double dt = 1.0;               // s  (time-step)
const double dx = 0.5;               // m  (segment length)

// Upstream boundary: constant concentration of 1 kg/m³ for the first 600 s,
// then zero (simulating a short-duration pulse).
int numTimeSteps = (int)(totalDuration / dt);
double[] upstreamBoundary = new double[numTimeSteps];
for (int i = 0; i < Math.Min(600, numTimeSteps); i++)
{
    upstreamBoundary[i] = 1.0; // kg/m³
}

// Physical / chemical parameters
var parameters = new StreamParameters(
    A: 1.0,        // main-channel cross-sectional area (m²)
    AS: 0.1,       // storage-zone cross-sectional area (m²)
    CL: 0.0,       // lateral inflow concentration (kg/m³)
    D: 0.05,       // dispersion coefficient (m²/s)
    Q: 0.5,        // flow rate (m³/s)
    QLin: 0.0,     // lateral inflow rate per unit length (m²/s)
    Alpha: 0.001,  // storage exchange rate (1/s)
    ChatS: 0.0,    // storage-zone equilibrium concentration (kg/m³)
    Kd: 0.0,       // partition coefficient (m³/kg)
    Lambda: 0.0,   // first-order decay (main channel) (1/s)
    LambdaS: 0.0,  // first-order decay (storage zone) (1/s)
    LambdaHat: 0.0,  // desorption rate from sediment (1/s)
    LambdaHatS: 0.0, // desorption rate in storage zone (1/s)
    Rho: 0.0);     // sediment bulk density (kg/m³)

var stream = new SoluteTransport.Stream(streamLength, totalDuration, dt, dx, upstreamBoundary, parameters);

Console.WriteLine("Running 1-D solute advection-dispersion simulation...");
Console.WriteLine($"  Stream length : {streamLength} m");
Console.WriteLine($"  Duration      : {totalDuration} s");
Console.WriteLine($"  Time step     : {dt} s");
Console.WriteLine($"  Segment length: {dx} m");
Console.WriteLine();

SimulationResult result = stream.SimulateAdvectionAndDispersion();

// Print the main-channel concentration at the mid-point and downstream end
// every 300 s.
int midSegment = (int)(streamLength / dx / 2);
int lastSegment = (int)(streamLength / dx) - 1;

Console.WriteLine($"{"Time (s)",10}  {"C mid (kg/m³)",15}  {"C end (kg/m³)",15}");
Console.WriteLine(new string('-', 46));

for (int j = 0; j < numTimeSteps; j += 300)
{
    Console.WriteLine($"{j * dt,10:F0}  {result.C[midSegment, j],15:F6}  {result.C[lastSegment, j],15:F6}");
}

Console.WriteLine();
Console.WriteLine("Simulation complete.");
