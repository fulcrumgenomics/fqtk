use anyhow::{Result, bail};
use enum_dispatch::enum_dispatch;

#[enum_dispatch]
pub trait Command {
    #[allow(clippy::missing_errors_doc)]
    fn execute(&self) -> Result<()>;
}

/// Returns an error listing every message in `errors`, which a command collects while checking
/// its arguments so that all the problems are reported at once, or `Ok` if there are none.
pub(crate) fn validation_result(errors: &[String]) -> Result<()> {
    if errors.is_empty() {
        return Ok(());
    }
    let details = errors.iter().fold(String::new(), |mut details, error| {
        details.push_str(&format!("    - {error}\n"));
        details
    });
    bail!("The following errors with the input(s) were detected:\n{details}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_errors_is_ok() {
        assert!(validation_result(&[]).is_ok());
    }

    #[test]
    fn every_error_is_listed() {
        let errors = ["First problem.".to_owned(), "Second problem.".to_owned()];
        let message = validation_result(&errors).unwrap_err().to_string();
        assert_eq!(
            message,
            "The following errors with the input(s) were detected:\n    - First problem.\n    - \
             Second problem.\n"
        );
    }
}
