# README  --  SNS-IK matlab code

This directory contains matlab implementations of a few of the basic SNS-IK
algorithms. These files should be considered experimental, and are primarily
included for reference and debugging purposes.

## Organization

This directory contains three types of files:
- `TEST_*_solvers.m` are entry-point scripts that are used to test either velocity or acceleration solvers.
- `runTest_*_.m` are test functions that operate on a specific class of solvers (*eg.* velocity or acceleration).
- `snsIk_*.m` and Matlab implementations for different variants of the SNS-IK solvers.
