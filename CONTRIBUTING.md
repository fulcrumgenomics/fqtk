
## Updating the README

The repository README contains a copy of the `fqtk demux` usage, which must be manually kept in sync with the tool.

When updating the tool's arguments or docstring, refresh the README by running the included script: 

```console
bash .github/scripts/update-docs.sh
```

GitHub Actions will verify that the README remains up-to-date, and commits with outdated usage will fail CI.

