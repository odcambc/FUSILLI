#!/usr/bin/env python3
"""
Unit tests for utils.py

Tests I/O functions and validation utilities.
"""

import pytest
import sys
from pathlib import Path

# Add workflow/scripts to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "workflow" / "scripts"))

from utils import (
    parse_fasta,
    parse_partners_csv,
    parse_samples_csv,
    validate_nucleotide_sequence,
    validate_sequences_match_config,
    write_fasta,
    write_csv,
)


# =============================================================================
# TEST: parse_fasta
# =============================================================================

class TestParseFasta:
    """Tests for FASTA file parsing."""

    def test_basic_parsing(self, tmp_path):
        """Should parse simple FASTA file."""
        fasta_content = """>seq1
ATGCATGC
>seq2
GCTAGCTA
"""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(fasta_content)

        sequences = parse_fasta(fasta_file)

        assert sequences['seq1'] == 'ATGCATGC'
        assert sequences['seq2'] == 'GCTAGCTA'

    def test_multiline_sequences(self, tmp_path):
        """Should handle multiline sequences."""
        fasta_content = """>seq1
ATGC
ATGC
ATGC
"""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(fasta_content)

        sequences = parse_fasta(fasta_file)

        assert sequences['seq1'] == 'ATGCATGCATGC'

    def test_uppercase_conversion(self, tmp_path):
        """Should convert sequences to uppercase."""
        fasta_content = """>seq1
atgcATGC
"""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(fasta_content)

        sequences = parse_fasta(fasta_file)

        assert sequences['seq1'] == 'ATGCATGC'

    def test_header_description_ignored(self, tmp_path):
        """Should use only first word of header as name."""
        fasta_content = """>seq1 this is a description
ATGC
"""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(fasta_content)

        sequences = parse_fasta(fasta_file)

        assert 'seq1' in sequences
        assert 'this' not in sequences

    def test_empty_lines_ignored(self, tmp_path):
        """Should skip empty lines."""
        fasta_content = """>seq1
ATGC

>seq2
GCTA
"""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(fasta_content)

        sequences = parse_fasta(fasta_file)

        assert len(sequences) == 2

    def test_file_not_found(self, tmp_path):
        """Should raise FileNotFoundError for missing file."""
        with pytest.raises(FileNotFoundError):
            parse_fasta(tmp_path / "nonexistent.fasta")

    def test_duplicate_name_error(self, tmp_path):
        """Should raise ValueError for duplicate sequence names."""
        fasta_content = """>seq1
ATGC
>seq1
GCTA
"""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(fasta_content)

        with pytest.raises(ValueError, match="Duplicate"):
            parse_fasta(fasta_file)

    def test_sequence_before_header_error(self, tmp_path):
        """Should raise ValueError for sequence before first header."""
        fasta_content = """ATGC
>seq1
GCTA
"""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(fasta_content)

        with pytest.raises(ValueError, match="before header"):
            parse_fasta(fasta_file)


# =============================================================================
# TEST: parse_partners_csv
# =============================================================================

class TestParsePartnersCsv:
    """Tests for partners CSV parsing."""

    def test_basic_parsing(self, tmp_path):
        """Should parse partners CSV file."""
        csv_content = """partner_name,sequence_length,include,description
TPR,426,true,TPR fusion partner
CCDC6,303,false,CCDC6 partner
"""
        csv_file = tmp_path / "partners.csv"
        csv_file.write_text(csv_content)

        partners = parse_partners_csv(csv_file)

        assert partners['TPR']['sequence_length'] == 426
        assert partners['TPR']['include'] == True
        assert partners['CCDC6']['include'] == False

    def test_include_variations(self, tmp_path):
        """Should handle various boolean representations."""
        csv_content = """partner_name,sequence_length,include,description
A,100,true,desc
B,100,True,desc
C,100,yes,desc
D,100,1,desc
E,100,false,desc
F,100,no,desc
"""
        csv_file = tmp_path / "partners.csv"
        csv_file.write_text(csv_content)

        partners = parse_partners_csv(csv_file)

        assert partners['A']['include'] == True
        assert partners['B']['include'] == True
        assert partners['C']['include'] == True
        assert partners['D']['include'] == True
        assert partners['E']['include'] == False
        assert partners['F']['include'] == False

    def test_skips_comment_lines(self, tmp_path):
        """Should skip lines starting with #."""
        csv_content = """# This is a comment
partner_name,sequence_length,include,description
TPR,426,true,desc
# Another comment
"""
        csv_file = tmp_path / "partners.csv"
        csv_file.write_text(csv_content)

        partners = parse_partners_csv(csv_file)

        assert 'TPR' in partners
        assert len(partners) == 1

    def test_missing_required_columns(self, tmp_path):
        """Should raise ValueError for missing required columns."""
        csv_content = """partner_name,include
TPR,true
"""
        csv_file = tmp_path / "partners.csv"
        csv_file.write_text(csv_content)

        with pytest.raises(ValueError, match="Missing required columns"):
            parse_partners_csv(csv_file)


# =============================================================================
# TEST: parse_samples_csv
# =============================================================================

class TestParseSamplesCsv:
    """Tests for samples CSV parsing."""

    def test_basic_parsing(self, tmp_path):
        """Should parse samples CSV file."""
        csv_content = """sample,condition,file
sample1,baseline,file1
sample2,treated,file2
"""
        csv_file = tmp_path / "samples.csv"
        csv_file.write_text(csv_content)

        samples = parse_samples_csv(csv_file)

        assert samples['sample1']['condition'] == 'baseline'
        assert samples['sample1']['file'] == 'file1'
        assert samples['sample2']['condition'] == 'treated'

    def test_optional_columns(self, tmp_path):
        """Should handle optional columns with defaults."""
        csv_content = """sample,condition,file,replicate,time,tile
sample1,baseline,file1,2,24h,A1
"""
        csv_file = tmp_path / "samples.csv"
        csv_file.write_text(csv_content)

        samples = parse_samples_csv(csv_file)

        assert samples['sample1']['replicate'] == 2
        assert samples['sample1']['time'] == '24h'
        assert samples['sample1']['tile'] == 'A1'

    def test_default_values(self, tmp_path):
        """Should use defaults for missing optional columns."""
        csv_content = """sample,condition,file
sample1,baseline,file1
"""
        csv_file = tmp_path / "samples.csv"
        csv_file.write_text(csv_content)

        samples = parse_samples_csv(csv_file)

        assert samples['sample1']['replicate'] == 1


# =============================================================================
# TEST: validate_nucleotide_sequence
# =============================================================================

class TestValidateNucleotideSequence:
    """Tests for nucleotide sequence validation."""

    def test_valid_sequence(self):
        """Should accept valid nucleotide sequences."""
        is_valid, error = validate_nucleotide_sequence("ATGCATGC")

        assert is_valid == True
        assert error is None

    def test_lowercase_valid(self):
        """Should accept lowercase (converted internally)."""
        is_valid, error = validate_nucleotide_sequence("atgc")

        assert is_valid == True

    def test_invalid_characters(self):
        """Should reject invalid characters."""
        is_valid, error = validate_nucleotide_sequence("ATGXYZ")

        assert is_valid == False
        assert 'X' in error

    def test_ambiguous_codes_rejected_by_default(self):
        """Should reject ambiguous codes by default."""
        is_valid, error = validate_nucleotide_sequence("ATGN")

        assert is_valid == False

    def test_ambiguous_codes_allowed(self):
        """Should accept ambiguous codes when allowed."""
        is_valid, error = validate_nucleotide_sequence("ATGN", allow_ambiguous=True)

        assert is_valid == True


# =============================================================================
# TEST: validate_sequences_match_config
# =============================================================================

class TestValidateSequencesMatchConfig:
    """Tests for configuration validation."""

    def test_valid_config(self):
        """Should return empty list for valid config."""
        sequences = {
            'Met_WT': 'ATGCATGC',
            'TPR': 'GCTAGCTA' * 50,  # 400 nt
        }
        partners = {
            'TPR': {'include': True, 'sequence_length': 400, 'description': ''},
        }

        messages = validate_sequences_match_config(sequences, partners, 'Met_WT')

        errors = [m for m in messages if m.startswith('ERROR')]
        assert len(errors) == 0

    def test_missing_anchor(self):
        """Should report error for missing anchor."""
        sequences = {
            'TPR': 'GCTAGCTA',
        }
        partners = {
            'TPR': {'include': True, 'sequence_length': 8, 'description': ''},
        }

        messages = validate_sequences_match_config(sequences, partners, 'Met_WT')

        assert any('Met_WT' in m and 'ERROR' in m for m in messages)

    def test_missing_partner(self):
        """Should report error for missing partner."""
        sequences = {
            'Met_WT': 'ATGCATGC',
        }
        partners = {
            'TPR': {'include': True, 'sequence_length': 8, 'description': ''},
        }

        messages = validate_sequences_match_config(sequences, partners, 'Met_WT')

        assert any('TPR' in m and 'ERROR' in m for m in messages)

    def test_length_mismatch_warning(self):
        """Should report warning for length mismatch."""
        sequences = {
            'Met_WT': 'ATGCATGC',
            'TPR': 'GCTAGCTA',  # 8 nt
        }
        partners = {
            'TPR': {'include': True, 'sequence_length': 100, 'description': ''},  # Wrong!
        }

        messages = validate_sequences_match_config(sequences, partners, 'Met_WT')

        assert any('mismatch' in m.lower() for m in messages)

    def test_excludes_non_included_partners(self):
        """Should skip validation for non-included partners."""
        sequences = {
            'Met_WT': 'ATGCATGC',
        }
        partners = {
            'TPR': {'include': False, 'sequence_length': 100, 'description': ''},
        }

        messages = validate_sequences_match_config(sequences, partners, 'Met_WT')

        # Should not error on TPR since it's not included
        assert not any('TPR' in m for m in messages)


# =============================================================================
# TEST: write_fasta
# =============================================================================

class TestWriteFasta:
    """Tests for FASTA file writing."""

    def test_basic_writing(self, tmp_path):
        """Should write valid FASTA file."""
        sequences = {
            'seq1': 'ATGCATGC',
            'seq2': 'GCTAGCTA',
        }
        output_file = tmp_path / "output.fasta"

        write_fasta(output_file, sequences)

        content = output_file.read_text()
        assert '>seq1' in content
        assert 'ATGCATGC' in content

    def test_creates_parent_dirs(self, tmp_path):
        """Should create parent directories if needed."""
        sequences = {'seq1': 'ATGC'}
        output_file = tmp_path / "subdir" / "output.fasta"

        write_fasta(output_file, sequences)

        assert output_file.exists()

    def test_line_wrapping(self, tmp_path):
        """Should wrap long sequences."""
        sequences = {'seq1': 'A' * 100}
        output_file = tmp_path / "output.fasta"

        write_fasta(output_file, sequences, line_width=50)

        lines = output_file.read_text().strip().split('\n')
        # Header + 2 sequence lines (100 / 50 = 2)
        assert len(lines) == 3


# =============================================================================
# TEST: write_csv
# =============================================================================

class TestWriteCsv:
    """Tests for CSV file writing."""

    def test_basic_writing(self, tmp_path):
        """Should write valid CSV file."""
        data = [
            {'name': 'A', 'value': 1},
            {'name': 'B', 'value': 2},
        ]
        output_file = tmp_path / "output.csv"

        write_csv(output_file, data, fieldnames=['name', 'value'])

        content = output_file.read_text()
        assert 'name,value' in content
        assert 'A,1' in content

    def test_header_comment(self, tmp_path):
        """Should add header comment if provided."""
        data = [{'name': 'A'}]
        output_file = tmp_path / "output.csv"

        write_csv(output_file, data, fieldnames=['name'],
                  header_comment="This is a comment")

        content = output_file.read_text()
        assert content.startswith('# This is a comment')


# =============================================================================
# RUN TESTS
# =============================================================================

if __name__ == '__main__':
    pytest.main([__file__, '-v'])
