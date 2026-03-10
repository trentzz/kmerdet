mod common;

use common::MockDb;
use kmerdet::walker::extension;

#[test]
fn test_extend_forward_filters_by_count() {
    let mut db = MockDb::new(4);
    // Set up a k-mer "ACGT" and its possible forward extensions
    db.set("ACGT", 100);
    db.set("CGTA", 50); // A extension: passes
    db.set("CGTC", 2);  // C extension: too low
    db.set("CGTG", 1);  // G extension: too low
    db.set("CGTT", 0);  // T extension: zero

    let children = extension::extend_forward(&db, "ACGT", 5, 0.05);

    assert_eq!(children.len(), 1);
    assert_eq!(children[0].sequence, "CGTA");
    assert_eq!(children[0].count, 50);
}

#[test]
fn test_extend_forward_all_pass() {
    let mut db = MockDb::new(4);
    db.set("CGTA", 100);
    db.set("CGTC", 100);
    db.set("CGTG", 100);
    db.set("CGTT", 100);

    let children = extension::extend_forward(&db, "ACGT", 2, 0.05);
    assert_eq!(children.len(), 4);
}

#[test]
fn test_extend_forward_none_pass() {
    let mut db = MockDb::new(4);
    // All zero counts
    let children = extension::extend_forward(&db, "ACGT", 2, 0.05);
    assert!(children.is_empty());
}
