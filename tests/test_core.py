import copy
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from io import BytesIO
from pathlib import Path
from unittest.mock import patch
import numpy as np
from pydantic import ValidationError
from omicsbench.models import Dataset, File, safe_relative
from omicsbench.registry import registry, lookup
from omicsbench.hashing import digest, contained
from omicsbench.cache import get, verify_cache
from omicsbench.download import transfer
from omicsbench.validate import validate
from omicsbench.rnaseq import ranked,correlation,check_design,contrast_selection_order,median_ratio_size_factors,select_genes,paired_sample,read_counts,preservation
from omicsbench.expression_atlas import load_design,load_design_fields,normalize_accession,parse_catalogue,stage

ROOT=Path(__file__).resolve().parents[1]

class CoreTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.root=Path(self.temp.name)
        shutil.copytree(ROOT/'datasets',self.root/'datasets')
        self.model,self.folder=lookup(self.root,'fixture-001')
    def tearDown(self):self.temp.cleanup()
    def test_manifest_round_trip(self):
        self.assertEqual(Dataset.model_validate_json(self.model.model_dump_json()),self.model)
    def test_required_families(self):
        for field in ['id','source','rights','reference','samples','derivation','files','validation','organism','taxon_id','release']:
            with self.subTest(field=field):
                data=self.model.model_dump(mode='json');del data[field]
                with self.assertRaises(ValidationError):Dataset.model_validate(data)
    def test_unknown_fields(self):
        data=self.model.model_dump();data['rights']['guessed']=True
        with self.assertRaises(ValidationError):Dataset.model_validate(data)
    def test_invalid_nested_values(self):
        for field,value in [('bytes',-1),('sha256','0'*63),('role','mystery'),('url','http://bad')]:
            data=self.model.files[0].model_dump();data[field]=value
            with self.subTest(field=field),self.assertRaises(ValidationError):File.model_validate(data)
    def test_duplicate_roles(self):
        data=self.model.model_dump();file=copy.deepcopy(data['files'][0]);file['path']='other.tsv';data['files'].append(file)
        # Use a data role to avoid intentionally repeatable documentation roles.
        data['files'][-1]['role']='raw_counts';data['files'][-1]['tier']='nano'
        with self.assertRaises(ValidationError):Dataset.model_validate(data)
    def test_duplicate_ids(self):
        shutil.copytree(self.folder,self.root/'datasets/duplicate/fixture-001')
        with self.assertRaisesRegex(ValueError,'Duplicate'):registry(self.root)
    def test_unknown_id(self):
        with self.assertRaisesRegex(ValueError,'fixture-001'):lookup(self.root,'fixture-002')
    def test_path_traversal(self):
        for path in ['../outside','/absolute','C:/file','a\\b','a/../b','./a','a//b']:
            with self.subTest(path=path),self.assertRaises(ValueError):safe_relative(path)
    def test_rights_gate(self):
        data=self.model.model_dump();data['rights']['status']='AMBER'
        with self.assertRaises(ValidationError):Dataset.model_validate(data)
    def test_validate(self):self.assertEqual(validate(self.model,self.folder)['status'],'pass')
    def test_corruption(self):
        path=self.folder/'nano/counts.tsv';body=path.read_bytes();path.write_bytes(body[:-1]+b'x')
        with self.assertRaisesRegex(ValueError,'checksum'):validate(self.model,self.folder)
    def test_undeclared(self):
        (self.folder/'hidden.tsv').write_text('hidden')
        with self.assertRaisesRegex(ValueError,'undeclared'):validate(self.model,self.folder)
    def test_sample_order(self):
        changed=self.model.model_copy(deep=True);changed.samples.reverse()
        with self.assertRaisesRegex(ValueError,'sample order'):validate(changed,self.folder)
    def test_confounded(self):
        samples=copy.deepcopy(self.model.samples)
        for sample in samples:sample.batch=sample.condition
        with self.assertRaisesRegex(ValueError,'confounded'):check_design(samples)
    def test_balanced_design(self):self.assertEqual(check_design(self.model.samples)[0].shape,(8,3))
    def test_known_ranks(self):self.assertEqual(ranked([3,1,1,4]).tolist(),[2,0.5,0.5,3])
    def test_known_correlations(self):
        self.assertAlmostEqual(correlation([1,2,3],[3,2,1]),-1)
        with self.assertRaises(ValueError):correlation([1,1,1],[1,2,3])
        with self.assertRaises(ValueError):correlation([1,2,np.nan],[1,2,3])
    def test_median_ratio_size_factors(self):
        counts=np.array([[10,20,40],[5,10,20],[0,3,7]])
        factors=median_ratio_size_factors(counts)
        self.assertTrue(np.allclose(factors,[0.5,1,2]))
        with self.assertRaises(ValueError):median_ratio_size_factors(np.array([[0,1],[1,0]]))
    def test_known_preservation_metrics(self):
        full={'effects':np.array([4.,3.,2.,1.]),'distances':np.array([1.,2.,3.])}
        pocket={'effects':np.array([-4.,-3.,-2.,-1.]),'distances':np.array([3.,2.,1.])}
        metrics=preservation(full,pocket,['a','b','c','d'],['a','b','c','d'],2)
        self.assertAlmostEqual(metrics['spearman_logfc'],-1)
        self.assertEqual(metrics['top_k_jaccard'],1)
        self.assertAlmostEqual(metrics['distance_correlation'],-1)
        self.assertEqual(metrics['sign_concordance'],0)
    def test_missing_file(self):
        (self.folder/'nano/counts.tsv').unlink()
        with self.assertRaisesRegex(ValueError,'missing'):validate(self.model,self.folder)
    def test_link_only_get(self):
        with self.assertRaisesRegex(ValueError,'link-only'):get(self.root,self.root/'cache','rnaseq-001','nano')
    def test_reference_gate(self):
        data=self.model.model_dump();data['reference']['compatibility']='unresolved'
        with self.assertRaises(ValidationError):Dataset.model_validate(data)
    def test_metric_threshold_enforced(self):
        m=self.model.model_copy(deep=True);m.validation.metrics[0].minimum=1.0
        with self.assertRaisesRegex(ValueError,'below'):validate(m,self.folder)
    def test_seeded_selection(self):
        genes=[str(i) for i in range(100)];counts=np.arange(800).reshape(100,8)
        self.assertEqual(select_genes(genes,counts,50,7),select_genes(genes,counts,50,7))
        self.assertNotEqual(select_genes(genes,counts,50,7),select_genes(genes,counts,50,8))
        self.assertEqual(select_genes(genes,counts,100,7),list(range(100)))
    def test_contrast_selection_anchors_effects(self):
        genes=[f'g{i:03}' for i in range(200)];counts=np.arange(1600).reshape(200,8);effects=np.arange(200,dtype=float)
        order=contrast_selection_order(genes,counts,effects,7,anchors=20)
        self.assertEqual(set(order[:20]),set(range(180,200)))
        self.assertEqual(order,contrast_selection_order(genes,counts,effects,7,anchors=20))
        self.assertEqual(len(set(order)),len(genes))
    def test_fractional_counts_rejected(self):
        p=self.root/'bad.tsv';p.write_text('gene_id\ta\tb\ng1\t1.5\t2\n')
        with self.assertRaisesRegex(ValueError,'integers'):read_counts(p)
    def test_cache_roundtrip(self):
        path=get(self.root,self.root/'cache','fixture-001','nano')
        self.assertTrue((path/'nano.complete.json').is_file())
        self.assertEqual(verify_cache(self.root,self.root/'cache')['verified_files'],1)
        (path/'nano/counts.tsv').write_text('bad')
        with self.assertRaisesRegex(ValueError,'corrupt'):verify_cache(self.root,self.root/'cache')
    def test_failed_download_atomic(self):
        target=self.root/'cache/file.tsv'
        with self.assertRaises(ValueError):transfer(self.folder/'nano/counts.tsv',target,'0'*64,10)
        self.assertFalse(target.exists());self.assertFalse(list(target.parent.glob('.partial-*')))
    def test_interrupted_download(self):
        target=self.root/'cache/file.tsv'
        class Broken:
            def geturl(self):return 'https://example.org/data'
            def __enter__(self):return self
            def __exit__(self,*args):pass
            def read(self,n):raise OSError('connection interrupted')
        with patch('urllib.request.urlopen',return_value=Broken()):
            with self.assertRaises(OSError):transfer('https://example.org/data',target,'0'*64,10)
        self.assertFalse(target.exists());self.assertFalse(list(target.parent.glob('.partial-*')))
    def test_failed_get_no_marker(self):
        (self.folder/'nano/counts.tsv').write_text('bad')
        with self.assertRaises(ValueError):get(self.root,self.root/'cache','fixture-001','nano')
        self.assertFalse(list((self.root/'cache').glob('**/*.complete.json')))
    def test_paired_reads(self):
        a=self.root/'r1.fastq';b=self.root/'r2.fastq'
        a.write_text('@x/1\nACT\n+\nIII\n@y/1\nAGT\n+\nIII\n')
        b.write_text('@x/2\nACT\n+\nIII\n@y/2\nAGT\n+\nIII\n')
        result=paired_sample(a,b,self.root/'o1',self.root/'o2',1,7,'sample')
        self.assertEqual(result['selected_pairs'],2)
        b.write_text('@wrong/2\nACT\n+\nIII\n')
        with self.assertRaisesRegex(ValueError,'read IDs'):paired_sample(a,b,self.root/'o1',self.root/'o2',1,7,'sample')
    def test_cli_commands(self):
        for command in [['list'],['info','fixture-001'],['validate','fixture-001'],['provenance','fixture-001'],['doctor'],['get','fixture-001'],['verify-cache']]:
            result=subprocess.run([sys.executable,'-m','omicsbench','--root',str(self.root),'--cache',str(self.root/'cache')]+command,capture_output=True,text=True)
            with self.subTest(command=command):
                self.assertEqual(result.returncode,0,result.stderr);json.loads(result.stdout)
    def test_cli_invalid_exit(self):
        result=subprocess.run([sys.executable,'-m','omicsbench','--root',str(self.root),'info','missing'],capture_output=True,text=True)
        self.assertEqual(result.returncode,1);self.assertNotIn('Traceback',result.stderr)
    def test_fixture_rebuild(self):
        destination=self.root/'rebuilt'
        result=subprocess.run([sys.executable,str(ROOT/'workflows/build_fixture.py'),str(destination)],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        rebuilt=destination/'datasets/rnaseq/fixture-001'
        for original in self.folder.rglob('*'):
            if original.is_file():self.assertEqual(digest(original),digest(rebuilt/original.relative_to(self.folder)))

    def test_atlas_accession(self):
        self.assertEqual(normalize_accession('e-mtab-8572'),'E-MTAB-8572')
        for value in ['E-MTAB-8572/other','8572','E-MTAB-x','E-../8572']:
            with self.subTest(value=value),self.assertRaises(ValueError):normalize_accession(value)

    def test_atlas_catalogue(self):
        rows=[
            {'type':'icon-raw-counts','description':'counts','url':'experiments-content/E-MTAB-8572/resources/counts'},
            {'type':'icon-experiment-design','description':'design','url':'experiments-content/E-MTAB-8572/resources/design'},
            {'type':'icon-tsv','description':'results','url':'experiments-content/E-MTAB-8572/resources/results'},
        ]
        found=parse_catalogue('E-MTAB-8572',json.dumps(rows).encode())
        self.assertEqual(set(found),{'raw_counts','experiment_design'})
        self.assertTrue(found['raw_counts'].url.startswith('https://www.ebi.ac.uk/gxa/'))

    def test_atlas_catalogue_rejects_cross_experiment(self):
        rows=[
            {'type':'icon-raw-counts','description':'counts','url':'experiments-content/E-MTAB-1/resources/counts'},
            {'type':'icon-experiment-design','description':'design','url':'experiments-content/E-MTAB-8572/resources/design'},
        ]
        with self.assertRaisesRegex(ValueError,'outside'):parse_catalogue('E-MTAB-8572',json.dumps(rows).encode())

    def test_atlas_stage_and_validate(self):
        accession='E-MTAB-8572'
        catalogue=json.dumps([
            {'type':'icon-raw-counts','description':'counts','url':f'experiments-content/{accession}/resources/counts'},
            {'type':'icon-experiment-design','description':'design','url':f'experiments-content/{accession}/resources/design'},
        ]).encode()
        counts=b'Gene ID\tGene Name\tRUN1\tRUN2\ng1\tA\t1\t2\ng2\tB\t0\t4\n'
        design=b'Run\tFactor Value[group]\tAnalysed\nRUN1\tcontrol\tYes\nRUN2\ttreated\tYes\n'
        class Response(BytesIO):
            def __init__(self,body,url):super().__init__(body);self.url=url
            def geturl(self):return self.url
            def __enter__(self):return self
            def __exit__(self,*args):self.close()
        def opener(request,timeout):
            url=request.full_url if hasattr(request,'full_url') else request
            body=catalogue if '/json/' in url else counts if url.endswith('/counts') else design
            return Response(body,url)
        result=stage(accession,self.root/'atlas',opener=opener)
        self.assertEqual((result['genes'],result['samples']),(2,2))
        self.assertEqual({r['role'] for r in result['resources']},{'raw_counts','experiment_design'})

    def test_atlas_stage_rejects_sample_mismatch(self):
        accession='E-MTAB-8572'
        catalogue=json.dumps([
            {'type':'icon-raw-counts','description':'counts','url':f'experiments-content/{accession}/resources/counts'},
            {'type':'icon-experiment-design','description':'design','url':f'experiments-content/{accession}/resources/design'},
        ]).encode()
        payloads={'counts':b'Gene ID\tGene Name\tRUN1\tRUN2\ng1\tA\t1\t2\n','design':b'Run\tAnalysed\nRUN1\tYes\nRUN3\tYes\n'}
        class Response(BytesIO):
            def __init__(self,body,url):super().__init__(body);self.url=url
            def geturl(self):return self.url
            def __enter__(self):return self
            def __exit__(self,*args):self.close()
        def opener(request,timeout):
            url=request.full_url if hasattr(request,'full_url') else request
            body=catalogue if '/json/' in url else payloads[url.rsplit('/',1)[-1]]
            return Response(body,url)
        with self.assertRaisesRegex(ValueError,'missing analysed'):stage(accession,self.root/'atlas-bad',opener=opener)

    def test_atlas_stage_records_unanalysed_count_columns(self):
        accession='E-MTAB-8572'
        catalogue=json.dumps([
            {'type':'icon-raw-counts','description':'counts','url':f'experiments-content/{accession}/resources/counts'},
            {'type':'icon-experiment-design','description':'design','url':f'experiments-content/{accession}/resources/design'},
        ]).encode()
        counts=b'Gene ID\tGene Name\tRUN1\tRUN2\tEXTRA\ng1\tA\t1\t2\t3\n'
        design=b'Run\tAnalysed\nRUN1\tYes\nRUN2\tYes\nEXTRA\tNo\n'
        class Response(BytesIO):
            def __init__(self,body,url):super().__init__(body);self.url=url
            def geturl(self):return self.url
            def __enter__(self):return self
            def __exit__(self,*args):self.close()
        def opener(request,timeout):
            url=request.full_url if hasattr(request,'full_url') else request
            body=catalogue if '/json/' in url else counts if url.endswith('/counts') else design
            return Response(body,url)
        result=stage(accession,self.root/'atlas-extra',opener=opener)
        self.assertEqual(result['samples'],2)
        self.assertEqual(result['raw_count_columns'],3)
        self.assertEqual(result['excluded_count_columns'],['EXTRA'])

    def test_atlas_design_fields(self):
        path=self.root/'design.tsv'
        path.write_text(
            'Run\tAnalysed\tFactor Value[group]\tSample Characteristic[individual]\n'
            'R1\tYes\tcontrol\tP1\nR2\tYes\ttreated\tP1\n'
            'R3\tYes\tcontrol\tP2\nR4\tYes\ttreated\tP2\n'
            'R5\tNo\tignored\tP3\n',encoding='utf-8'
        )
        fields=load_design_fields(path,['Factor Value[group]','Sample Characteristic[individual]'])
        self.assertEqual(fields['R2'],{'Factor Value[group]':'treated','Sample Characteristic[individual]':'P1'})
        self.assertNotIn('R5',fields)
        self.assertEqual(load_design(path,'Factor Value[group]')['R3'],'control')

    def test_atlas_design_rejects_duplicate_runs(self):
        path=self.root/'duplicate-design.tsv'
        path.write_text(
            'Run\tAnalysed\tFactor Value[group]\n'
            'R1\tYes\tcontrol\nR1\tYes\ttreated\nR2\tYes\tcontrol\nR3\tYes\ttreated\n',
            encoding='utf-8'
        )
        with self.assertRaisesRegex(ValueError,'unique'):
            load_design_fields(path,['Factor Value[group]'])

if __name__=='__main__':unittest.main()
