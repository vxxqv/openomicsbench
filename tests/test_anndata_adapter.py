import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch

from omicsbench.anndata_adapter import export_h5ad


class AnnDataAdapterTests(unittest.TestCase):
    def test_selected_columns_and_embedding_are_exported(self):
        example = types.SimpleNamespace(obs_names=['cell-A', 'cell-B'],
                                        obs={'kind': ['T', 'B'], 'batch': ['one', 'two']},
                                        obsm={'X_umap': [[0.0, 1.0], [2.0, 3.0]]}, file=None)
        module = types.SimpleNamespace(read_h5ad=lambda source, backed: example)
        with tempfile.TemporaryDirectory() as temporary, patch.dict(sys.modules, {'anndata': module}):
            directory = Path(temporary) / 'export'
            report = export_h5ad(Path(temporary) / 'example.h5ad', directory,
                                 label_column='kind', batch_column='batch', embedding_key='X_umap')
            self.assertEqual(report['cells'], 2)
            self.assertIn('cell-A\tT\tone', (directory / 'labels.tsv').read_text(encoding='utf-8'))
            self.assertIn('cell-B\t2.0\t3.0', (directory / 'embedding.tsv').read_text(encoding='utf-8'))
            with self.assertRaisesRegex(ValueError, 'already exists'):
                export_h5ad(Path(temporary) / 'example.h5ad', directory, label_column='kind')


if __name__ == '__main__':
    unittest.main()
