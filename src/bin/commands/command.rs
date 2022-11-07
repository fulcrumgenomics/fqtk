use anyhow::Result;
use enum_dispatch::enum_dispatch;

#[enum_dispatch]
pub trait Command {
    #[allow(clippy::missing_errors_doc)]
    fn execute(&self) -> Result<()>;
}
