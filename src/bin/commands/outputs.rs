//! Output files for the subcommands that write BGZF FASTQs: checking output paths before
//! anything is opened, creating pooled BGZF writers, and removing incomplete outputs after a
//! failure.

use anyhow::{Context, Result};
use log::warn;
use pooled_writer::{Pool, PoolBuilder, PooledWriter, bgzf::BgzfCompressor};
use std::fs::{self, File};
use std::io::BufWriter;
use std::path::{Path, PathBuf};

/// Checks `outputs` before any of them is created, adding a message to `errors` for an output
/// whose directory doesn't exist, that is also one of `inputs`, or that repeats an earlier
/// output.
pub(crate) fn check_output_paths(
    inputs: &[PathBuf],
    outputs: &[PathBuf],
    errors: &mut Vec<String>,
) {
    // Creating an output truncates it, so one that is also an input would be emptied before
    // it is read.  Comparing resolved paths catches different spellings of the same file.
    let resolved_inputs: Vec<PathBuf> =
        inputs.iter().filter_map(|input| resolved_path(input)).collect();
    let mut resolved_outputs: Vec<PathBuf> = Vec::with_capacity(outputs.len());
    for output in outputs {
        if let Some(parent) = output.parent().filter(|p| !p.as_os_str().is_empty() && !p.exists()) {
            // Outputs often share a directory, which only needs reporting once.
            let message = format!("Output parent directory {parent:?} does not exist.");
            if !errors.contains(&message) {
                errors.push(message);
            }
            continue;
        }
        let Some(resolved) = resolved_path(output) else { continue };
        if resolved_inputs.contains(&resolved) {
            errors.push(format!("Output {output:?} is also an input."));
        } else if resolved_outputs.contains(&resolved) {
            errors.push(format!("Output {output:?} is given more than once."));
        }
        resolved_outputs.push(resolved);
    }
}

/// Creates a BGZF writer for each of `paths`, in order, sharing a pool of `threads` compression
/// threads that compress at `compression_level`.  Returns the pool and the writers.
pub(crate) fn create_bgzf_writers(
    paths: &[PathBuf],
    threads: usize,
    compression_level: u8,
) -> Result<(Pool, Vec<PooledWriter>)> {
    let mut pool_builder = PoolBuilder::<_, BgzfCompressor>::new()
        .threads(threads)
        .queue_size(threads * 50)
        .compression_level(compression_level)?;
    let mut writers = Vec::with_capacity(paths.len());
    for path in paths {
        let file =
            File::create(path).with_context(|| format!("Failed to create output file {path:?}"))?;
        writers.push(pool_builder.exchange(BufWriter::new(file)));
    }
    Ok((pool_builder.build()?, writers))
}

/// Deletes each of `paths` that is a regular file, so that a failed run doesn't leave output
/// that looks complete.  Pipes, devices and symlinks are left alone: deleting a symlink removes
/// the link rather than its target, and on Linux `/dev/stdout` is one.
pub(crate) fn remove_incomplete_outputs(paths: &[PathBuf]) {
    for path in paths {
        if fs::symlink_metadata(path).is_ok_and(|m| m.is_file()) && fs::remove_file(path).is_ok() {
            warn!("Removed incomplete output {path:?}.");
        }
    }
}

/// Returns `path` with symlinks and relative components resolved, whether or not the file
/// exists yet, or `None` if its directory can't be resolved.
fn resolved_path(path: &Path) -> Option<PathBuf> {
    if path.exists() {
        return fs::canonicalize(path).ok();
    }
    let parent = path.parent().filter(|p| !p.as_os_str().is_empty()).unwrap_or(Path::new("."));
    Some(fs::canonicalize(parent).ok()?.join(path.file_name()?))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Read, Write};
    use tempfile::TempDir;

    /// Returns the messages `check_output_paths` adds for `inputs` and `outputs`.
    fn output_errors(inputs: &[PathBuf], outputs: &[PathBuf]) -> Vec<String> {
        let mut errors = Vec::new();
        check_output_paths(inputs, outputs, &mut errors);
        errors
    }

    #[test]
    fn new_outputs_in_an_existing_directory_are_accepted() {
        let dir = TempDir::new().unwrap();
        let input = dir.path().join("in.fq");
        fs::write(&input, "").unwrap();

        let outputs = [dir.path().join("out.R1.fq.gz"), dir.path().join("out.R2.fq.gz")];
        assert!(output_errors(&[input], &outputs).is_empty());
    }

    #[test]
    fn an_output_that_is_an_input_is_rejected() {
        let dir = TempDir::new().unwrap();
        let input = dir.path().join("in.fq.gz");
        fs::write(&input, "").unwrap();

        let errors = output_errors(&[input.clone()], &[input]);
        assert_eq!(errors.len(), 1);
        assert!(errors[0].contains("is also an input"), "{errors:?}");
    }

    #[test]
    fn an_output_that_is_an_input_spelled_differently_is_rejected() {
        let dir = TempDir::new().unwrap();
        fs::create_dir(dir.path().join("sub")).unwrap();
        let input = dir.path().join("in.fq.gz");
        fs::write(&input, "").unwrap();

        let output = dir.path().join("sub").join("..").join("in.fq.gz");
        assert!(output_errors(&[input], &[output])[0].contains("is also an input"));
    }

    #[test]
    fn an_output_given_twice_is_rejected() {
        let dir = TempDir::new().unwrap();
        let output = dir.path().join("out.fq.gz");

        let errors = output_errors(&[], &[output.clone(), output]);
        assert_eq!(errors.len(), 1);
        assert!(errors[0].contains("is given more than once"), "{errors:?}");
    }

    #[test]
    fn a_missing_output_directory_is_reported_once_for_all_its_outputs() {
        let dir = TempDir::new().unwrap();
        let missing = dir.path().join("missing");
        let outputs = [missing.join("out.R1.fq.gz"), missing.join("out.R2.fq.gz")];

        let errors = output_errors(&[], &outputs);
        assert_eq!(errors.len(), 1);
        assert!(errors[0].contains("does not exist"), "{errors:?}");
    }

    #[test]
    fn writers_write_bgzf_to_their_paths_in_order() {
        let dir = TempDir::new().unwrap();
        let paths = [dir.path().join("a.fq.gz"), dir.path().join("b.fq.gz")];
        let (mut pool, writers) = create_bgzf_writers(&paths, 2, 1).unwrap();
        for (mut writer, text) in writers.into_iter().zip(["first\n", "second\n"]) {
            writer.write_all(text.as_bytes()).unwrap();
            writer.close().unwrap();
        }
        pool.stop_pool().unwrap();

        for (path, expected) in paths.iter().zip(["first\n", "second\n"]) {
            let bytes = fs::read(path).unwrap();
            assert_eq!(&bytes[12..14], b"BC", "{path:?} is not BGZF");
            let mut text = String::new();
            flate2::read::MultiGzDecoder::new(&bytes[..]).read_to_string(&mut text).unwrap();
            assert_eq!(text, expected);
        }
    }

    #[test]
    fn an_output_that_cannot_be_created_is_an_error_naming_it() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("missing").join("out.fq.gz");

        let err = create_bgzf_writers(&[path.clone()], 2, 1).err().unwrap();
        assert!(format!("{err:#}").contains(&format!("{path:?}")), "{err:#}");
    }

    #[test]
    fn regular_file_outputs_are_removed() {
        let dir = TempDir::new().unwrap();
        let output = dir.path().join("out.fq.gz");
        fs::write(&output, "partial").unwrap();

        remove_incomplete_outputs(&[output.clone()]);
        assert!(!output.exists());
    }

    #[test]
    fn missing_outputs_are_ignored() {
        let dir = TempDir::new().unwrap();
        remove_incomplete_outputs(&[dir.path().join("never-created.fq.gz")]);
    }

    #[cfg(unix)]
    #[test]
    fn a_symlinked_output_and_its_target_are_left_in_place() {
        let dir = TempDir::new().unwrap();
        let target = dir.path().join("target.fq.gz");
        fs::write(&target, "partial").unwrap();
        let link = dir.path().join("out.fq.gz");
        std::os::unix::fs::symlink(&target, &link).unwrap();

        remove_incomplete_outputs(&[link.clone()]);
        assert!(fs::symlink_metadata(&link).unwrap().file_type().is_symlink());
        assert!(target.exists());
    }
}
