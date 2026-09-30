"""Gut mycobiome: fungal sequence signal, supported fungi, strains, and the
experimental Mycobiome Health Score (MHS-E1).

Implements BUILD_SPEC_v08.2. The module is additive: nothing here changes a
bacterial number, and a fungal percentage is never summed into a bacterial
composition. Every quantity carries its denominator; every score
contribution carries its source, cap and activation; every unresolved
strain question carries its reason.
"""
