// Path quantification (NNLS) tests.
// TODO Phase 4: Add tests when quantifier is implemented.

#[test]
#[ignore = "quantifier not yet implemented"]
fn test_quantify_single_path() {
    // Trivial case: one path should get rVAF = 1.0
}

#[test]
#[ignore = "quantifier not yet implemented"]
fn test_quantify_two_paths_equal() {
    // Two paths with equal counts should get rVAF ≈ 0.5 each
}

#[test]
#[ignore = "quantifier not yet implemented"]
fn test_quantify_negative_coefficients() {
    // NNLS should clip negative coefficients to zero
}
