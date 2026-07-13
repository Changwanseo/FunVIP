# FunVIP tests

Fast, dependency-light unit tests for the pure/near-pure functions, plus (marked)
slower integration checks.

```bash
pip install -e ".[test]"      # or: conda run -n FunVIP_dev pip install pytest
pytest                        # unit tests only (no external tools needed)
pytest -m integration         # end-to-end (needs the bundled external tools + --email)
```

The unit tests (`test_hasher`, `test_logics`, `test_tool`, `test_validate_option`,
`test_cluster`) require none of the external aligners/tree tools and run in
seconds, so they are safe for CI. Integration tests that shell out to
BLAST/mmseqs/MAFFT/FastTree are marked `@pytest.mark.integration` and skipped by
default.
