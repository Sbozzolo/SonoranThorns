# SonoranThorns

`SonoranThorns` is a collection of thorns for the [Einstein
Toolkit](https://einsteintoolkit.org/).

## Christoffel

`Christoffel` computes the four-dimensional Christoffel symbols using the
derivatives of the metric. We compute the spatial derivatives with finite
difference methods (4th or 6th order), and we read off the time derivatives from
the right-hand-sides of the evolution equations.

`Christoffel` was tested by checking the output with the expected values
computed analytically for a boosted Kerr-Schild black hole. It agrees well and
converges at the expected rate.

The test in `Christoffel` requires a patched version of `LeanBSSNMoL` that
computes the right-hand-side of the evolution equations after the initial data.

`Christoffel` optionally saves and compute the derivatives of the
four-dimensional metric.

`Christoffel` defined up to 80 grid functions with three timelevels, so it can
be very memory-intensive!
