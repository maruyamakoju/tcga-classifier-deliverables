## Summary

-

## Type of change

- [ ] Code behavior
- [ ] Documentation
- [ ] Release artifact or metadata
- [ ] Tests / CI

## Checks

- [ ] `python -m pytest -q -rs`
- [ ] `python audit_publication_readiness.py`
- [ ] `python run_release_acceptance.py --timeout-seconds 300`
- [ ] If release files changed: `python build_release_lite.py --smoke --timeout-seconds 300`

## Safety

- [ ] No credentials, PHI, private sample IDs, or unpublished raw data are included.
- [ ] The research-only and GDC STAR-Counts input boundary remains clear.
