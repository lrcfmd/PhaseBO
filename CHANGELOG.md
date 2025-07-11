# Changelog

All notable changes to this project will be documented in this file.

---

## [0.1.4] - 2025-07-09
### Changed
- Migrated from `setup.py` to `pyproject.toml` for modern Python packaging.
- Updated to follow standard Semantic Versioning (three-part version).
- Improved documentation and clarified usage examples for crystal structure exploration workflows.

### Fixed
- Minor packaging and dependency fixes for smoother local editable installs.

---

## [0.1.3.1] - 2022-07-23
### Added
- Preliminary support for compositional space exploration with explicit DFT and interatomic potential integration.
- Early support for parallel or batch candidate evaluations.

### Changed
- Refined Thompson Sampling and Expected Improvement (EI) acquisition function implementations for more robust convergence.

---

## [0.1.3] - 2022-03-10
### Added
- Initial integration with external crystal structure prediction modules
- New plotting utilities to visualise acquisition surfaces and posterior distributions.

---

## [0.1.2] - 2021-10-07
### Added
- Support for multiple acquisition strategies: Expected Improvement (EI) and Thompson Sampling.

### Changed
- Improved Gaussian Process model initialisation for better uncertainty handling in high-dimensional compositional spaces.

---

## [0.1.1] - 2021-08-21
### Added
- Basic Bayesian optimisation loop for accelerated compositional exploration.
- Early demonstration examples on phase fields.

---

## [0.1.0] - 2021-03-05
### Added
- Initial prototype release.
- Core Bayesian optimisation engine with Gaussian Process support.
- Simple compositional candidate generator and initial evaluation interface.

---

## [Unreleased]
### Planned
- Support for multi-objective acquisition strategies.
- Integration with active learning for iterative interatomic potential refinement.
- Additional crystal structure filtering and ranking tools.
