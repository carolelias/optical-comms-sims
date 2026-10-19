# Optical communications simulations (MATLAB)

Four simulations written for the course *SEL0366 – Optical Communications* at the
University of São Paulo (EESC-USP). Each script is self-contained: it sets its
parameters, runs the model, plots the results and ends with a commented analysis.

> Code comments and analyses are in Portuguese, the language the course was taught in.

| Script | Topic | Method |
| --- | --- | --- |
| `FiltroDBR_CarolinaElias.m` | Narrow-band DBR filter at 1550 nm (TiO₂/SiO₂ Bragg mirrors with a half-wave defect) | Transfer-matrix method, TE/TM at 30° incidence |
| `ARROW_CarolinaElias.m` | Leakage loss of a waveguide over an antiresonant Si mirror | Leaky-mode solver (complex β) built on the TMM |
| `SplitStep_CarolinaElias.m` | Soliton propagation and third-order dispersion in silica fibre | Symmetric split-step Fourier method (NLSE) |
| `MonerieFibraW_CarolinaElias.m` | Dispersion-flattened W fibre (Monerie, 1982) | LP₀₁ dispersion equation (4×4 Bessel determinant) + Gloge waveguide dispersion |

## Summary

**1. DBR filter.** Spectral response of an (HL)ᴺ H₂ (LH)ᴺ microcavity for TE and TM
polarisation, plotted in wavelength (slightly asymmetric) and frequency (symmetric,
since the phase is linear in ω for dispersionless media). The Q factor grows roughly
exponentially with the number of periods (Q_TE ≈ 1.5·10³ at N = 5, ≈ 1.5·10⁶ at N = 11)
and is always higher for TE, whose interface reflectivity increases at oblique incidence.
Half-maximum points are found by bisection, so arbitrarily narrow resonances are resolved.

**2. ARROW.** A high-index core (n = 4.0) is leaky towards a Si substrate through a finite
cladding. Inserting a Si layer of thickness *t* in the cladding gives sharp loss peaks
(≈ 1000× the bare-guide loss) at the resonances k_x t = mπ + 2 arctan(γ/k_x), i.e. when the
mode phase-matches a guided mode of the Si layer, while at antiresonance the layer is
optically invisible (loss ≈ 0.98× the bare guide). The leaky mode is found as a complex
root of the TMM dispersion function, with mode tracking through anticrossings by field
confinement. A lossless and an absorbing mirror (n = 3.6 − 0.2j) are compared.

**3. Split-step Fourier.** β₂ and β₃ of silica at 1550 nm are obtained from the Sellmeier
equation (β₂ = −27.9 ps²/km, D = 21.9 ps/(nm·km), zero-dispersion wavelength 1273 nm),
including a check of how the polynomial order biases β₃ when fitting β(ω). Three scenarios
in soliton units:
- N = 1 fundamental soliton over 10 L_D: shape preserved to ~10⁻⁴, energy conserved (validation);
- N = 3 higher-order soliton: periodic evolution with period π/2;
- N = 2 with the real β₃ of silica (δ₃ ≈ 0.018, 50 fs pulse): soliton fission and emission
  of a blue-shifted dispersive wave at (ω − ω₀)T₀ ≈ +28, close to the phase-matching
  estimate 1/(2δ₃). This scenario uses a ±1024 T₀ window so the dispersive wave does not
  wrap around the periodic FFT boundary.

Raman scattering and self-steepening are not included.

**4. W fibre.** Universal B(V) curves of the LP₀₁ mode for several index-step ratios R,
showing the finite cutoff of the fundamental mode in depressed-cladding fibres. For
profile (b) of Monerie's Fig. 9 (2a = 6.5 µm, b/a = 1.5, R = −0.5), material and
waveguide dispersion cancel around 1.5–1.55 µm: the total dispersion stays within
~1 ps/(nm·km) between ~1.46 and 1.60 µm.

## Running

Open a script in MATLAB and run it. Local functions inside scripts require **R2016b or
later**; `xline`, `yline` and `sgtitle` require **R2018b or later**. No toolboxes are needed.
The ARROW sweep and the wide-window split-step scenario take from a few seconds to a
minute.

## References

- G. P. Agrawal, *Nonlinear Fiber Optics*, 5th ed., Academic Press (2013) — NLSE and split-step method (Sec. 2.4.1, App. B).
- M. A. Duguay, Y. Kokubun, T. L. Koch, L. Pfeiffer, "Antiresonant reflecting optical waveguides in SiO₂–Si multilayer structures", *Appl. Phys. Lett.* 49, 13 (1986).
- M. Monerie, "Propagation in doubly clad single-mode fibers", *IEEE J. Quantum Electron.* 18, 535 (1982).
- D. Gloge, "Dispersion in weakly guiding fibers", *Appl. Opt.* 10, 2442 (1971).
- S. J. Orfanidis, *Electromagnetic Waves and Antennas* — transfer matrices for multilayer structures.

Author: Carolina Elias de Almeida Américo
