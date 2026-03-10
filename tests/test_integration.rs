use assert_cmd::Command;

#[test]
fn test_cli_help() {
    Command::cargo_bin("kmerdet")
        .unwrap()
        .arg("--help")
        .assert()
        .success()
        .stdout(predicates::str::contains("kmerdet detects variants"))
        .stdout(predicates::str::contains("detect"))
        .stdout(predicates::str::contains("filter"))
        .stdout(predicates::str::contains("merge"))
        .stdout(predicates::str::contains("stats"))
        .stdout(predicates::str::contains("plot"))
        .stdout(predicates::str::contains("coverage"))
        .stdout(predicates::str::contains("run"));
}

#[test]
fn test_cli_version() {
    Command::cargo_bin("kmerdet")
        .unwrap()
        .arg("--version")
        .assert()
        .success()
        .stdout(predicates::str::contains("kmerdet"));
}

#[test]
fn test_detect_help() {
    Command::cargo_bin("kmerdet")
        .unwrap()
        .args(["detect", "--help"])
        .assert()
        .success()
        .stdout(predicates::str::contains("Jellyfish database"))
        .stdout(predicates::str::contains("--targets"));
}

#[test]
fn test_filter_help() {
    Command::cargo_bin("kmerdet")
        .unwrap()
        .args(["filter", "--help"])
        .assert()
        .success()
        .stdout(predicates::str::contains("--use-alt"))
        .stdout(predicates::str::contains("--min-coverage"));
}

#[test]
fn test_run_help() {
    Command::cargo_bin("kmerdet")
        .unwrap()
        .args(["run", "--help"])
        .assert()
        .success()
        .stdout(predicates::str::contains("Full pipeline"));
}
