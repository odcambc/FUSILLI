"""
FUSILLI Test Suite

Tests for the FUSILLI (Fusion Utility for Scanning and Identification of
Library Linked Interactions) pipeline.

Test Organization:
- test_fusion_sequences.py: Unit tests for breakpoint generation
- test_string_matcher.py: Unit tests for fusion detection
- test_utils.py: Unit tests for I/O and validation utilities
- test_integration.py: End-to-end integration tests

Run tests with:
    pytest                      # Run all tests
    pytest -v                   # Verbose output
    pytest tests/test_fusion_sequences.py  # Run specific file
    pytest -m "not slow"       # Skip slow tests
    pytest -m integration      # Run only integration tests
"""
