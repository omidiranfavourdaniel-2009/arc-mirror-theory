#!/usr/bin/env python3
"""First-order thermal-expansion sensitivity reference.

This is deliberately labelled a sensitivity model, not physical validation.
The final 2 m optical baseline is 0.01215 mrad; the 0.037 mrad value belongs
to the separate 200 mm reference case and is not used here.
"""

D_M=2.0
CTE=23.6e-6
BASELINE_SLOPE_RAD=0.01215e-3

print("Arc-Mirror thermal sensitivity reference")
print(f"Diameter: {D_M:.1f} m")
print(f"Aluminium free-expansion scale per kelvin: {CTE*D_M*1e6:.2f} micrometres/K")
print(f"Optical baseline slope error: {BASELINE_SLOPE_RAD*1e3:.5f} mrad")
print("Model status: first-order sensitivity analysis; not experimental validation.")
print("Dominant physical validation risks: thermal drift, thin-shell flex, machining accuracy and gravity/wind deflection.")
