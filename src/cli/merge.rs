use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;

use anyhow::{Context, Result};

#[derive(clap::Args, Debug)]
pub struct MergeArgs {
    /// Input result files to merge
    #[arg(short, long, num_args = 1..)]
    pub input: Vec<PathBuf>,

    /// Deduplicate rows with identical key columns
    #[arg(long)]
    pub deduplicate: bool,
}

pub fn run(args: MergeArgs, global: &super::GlobalOptions) -> Result<()> {
    anyhow::ensure!(!args.input.is_empty(), "at least one input file required");

    let mut output: Box<dyn Write> = match &global.output {
        Some(path) => Box::new(
            File::create(path)
                .with_context(|| format!("creating output file: {}", path.display()))?,
        ),
        None => Box::new(std::io::stdout().lock()),
    };

    let mut expected_header: Option<String> = None;
    let mut seen_keys: HashSet<String> = HashSet::new();

    for path in &args.input {
        let file = File::open(path)
            .with_context(|| format!("opening input file: {}", path.display()))?;
        let reader = BufReader::new(file);

        for (i, line) in reader.lines().enumerate() {
            let line = line?;
            if i == 0 {
                // Header line
                match &expected_header {
                    None => {
                        expected_header = Some(line.clone());
                        if !global.no_header {
                            writeln!(output, "{}", line)?;
                        }
                    }
                    Some(expected) => {
                        if &line != expected {
                            anyhow::bail!(
                                "header mismatch in {}: expected '{}', got '{}'",
                                path.display(),
                                expected,
                                line
                            );
                        }
                        // Skip duplicate header from subsequent files
                    }
                }
            } else {
                // Data line
                if args.deduplicate {
                    // Use sample + target + variant_name as dedup key
                    // Fields: sample(0), target(1), type(2), variant_name(3), ...
                    let fields: Vec<&str> = line.split('\t').collect();
                    let key = format!(
                        "{}|{}|{}",
                        fields.first().unwrap_or(&""),
                        fields.get(1).unwrap_or(&""),
                        fields.get(3).unwrap_or(&"")
                    );
                    if !seen_keys.insert(key) {
                        continue; // Duplicate, skip
                    }
                }
                writeln!(output, "{}", line)?;
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn write_test_tsv(dir: &std::path::Path, name: &str, content: &str) -> PathBuf {
        let path = dir.join(name);
        let mut f = File::create(&path).unwrap();
        f.write_all(content.as_bytes()).unwrap();
        path
    }

    #[test]
    fn test_merge_two_files() {
        let dir = tempfile::tempdir().unwrap();
        let header = "sample\ttarget\ttype\tvariant_name\trVAF";
        let f1 = write_test_tsv(
            dir.path(),
            "a.tsv",
            &format!("{}\ns1\tTP53\tSubstitution\tv1\t0.05\n", header),
        );
        let f2 = write_test_tsv(
            dir.path(),
            "b.tsv",
            &format!("{}\ns2\tBRCA1\tDeletion\tv2\t0.10\n", header),
        );

        let out_path = dir.path().join("merged.tsv");
        let global = super::super::GlobalOptions {
            config: None,
            threads: 0,
            verbose: 0,
            quiet: false,
            output: Some(out_path.clone()),
            format: super::super::OutputFormat::Tsv,
            no_header: false,
        };
        let args = MergeArgs {
            input: vec![f1, f2],
            deduplicate: false,
        };

        run(args, &global).unwrap();

        let content = std::fs::read_to_string(&out_path).unwrap();
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines.len(), 3); // header + 2 data lines
        assert!(lines[0].starts_with("sample\t"));
        assert!(lines[1].contains("s1"));
        assert!(lines[2].contains("s2"));
    }

    #[test]
    fn test_merge_deduplicate() {
        let dir = tempfile::tempdir().unwrap();
        let header = "sample\ttarget\ttype\tvariant_name\trVAF";
        let f1 = write_test_tsv(
            dir.path(),
            "a.tsv",
            &format!("{}\ns1\tTP53\tSubstitution\tv1\t0.05\n", header),
        );
        let f2 = write_test_tsv(
            dir.path(),
            "b.tsv",
            &format!("{}\ns1\tTP53\tSubstitution\tv1\t0.07\n", header),
        );

        let out_path = dir.path().join("merged.tsv");
        let global = super::super::GlobalOptions {
            config: None,
            threads: 0,
            verbose: 0,
            quiet: false,
            output: Some(out_path.clone()),
            format: super::super::OutputFormat::Tsv,
            no_header: false,
        };
        let args = MergeArgs {
            input: vec![f1, f2],
            deduplicate: true,
        };

        run(args, &global).unwrap();

        let content = std::fs::read_to_string(&out_path).unwrap();
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines.len(), 2); // header + 1 data line (duplicate removed)
        assert!(lines[1].contains("0.05")); // first occurrence wins
    }

    #[test]
    fn test_merge_header_mismatch() {
        let dir = tempfile::tempdir().unwrap();
        let f1 = write_test_tsv(dir.path(), "a.tsv", "col_a\tcol_b\nv1\tv2\n");
        let f2 = write_test_tsv(dir.path(), "b.tsv", "col_x\tcol_y\nv3\tv4\n");

        let global = super::super::GlobalOptions {
            config: None,
            threads: 0,
            verbose: 0,
            quiet: false,
            output: None,
            format: super::super::OutputFormat::Tsv,
            no_header: false,
        };
        let args = MergeArgs {
            input: vec![f1, f2],
            deduplicate: false,
        };

        let result = run(args, &global);
        assert!(result.is_err());
        assert!(
            result
                .unwrap_err()
                .to_string()
                .contains("header mismatch")
        );
    }

    #[test]
    fn test_merge_empty_input() {
        let global = super::super::GlobalOptions {
            config: None,
            threads: 0,
            verbose: 0,
            quiet: false,
            output: None,
            format: super::super::OutputFormat::Tsv,
            no_header: false,
        };
        let args = MergeArgs {
            input: vec![],
            deduplicate: false,
        };

        let result = run(args, &global);
        assert!(result.is_err());
        assert!(
            result
                .unwrap_err()
                .to_string()
                .contains("at least one input file")
        );
    }
}
