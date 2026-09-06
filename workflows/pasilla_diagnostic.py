"""Stage a pinned Pasilla source and evaluate reductions locally."""
import argparse,csv,json,platform
from pathlib import Path
import numpy as np
from omicsbench.download import transfer
from omicsbench.hashing import digest
from omicsbench.registry import lookup
from omicsbench.rnaseq import read_counts,diagnostic,preservation,select_genes
from build_fixture import write_counts,write_json

ROOT=Path(__file__).resolve().parents[1]
def run(staging):
    config=json.loads((ROOT/'workflows/pasilla-source.json').read_text())
    model,_=lookup(ROOT,config['dataset_id'])
    source=staging/'source';source.mkdir(parents=True,exist_ok=True)
    for f in config['sources']:
        transfer(f['url'],source/Path(f['file']).name,f['sha256'],f['bytes'])
    genes,ids,counts=read_counts(source/'pasilla_gene_counts.tsv')
    if ids != [s.sample_id for s in model.samples]:raise ValueError('Source sample order changed')
    full=diagnostic(counts,model.samples)
    source_size=(source/'pasilla_gene_counts.tsv').stat().st_size
    curve=[];selected=None
    for size in config['candidate_sizes']:
        indices=select_genes(genes,counts,size,config['seed'])
        names=[genes[i] for i in indices];pocket=diagnostic(counts[indices],model.samples)
        metrics=preservation(full,pocket,genes,names,config['top_k'])
        path=staging/f'candidates/{size}/counts.tsv'
        write_counts(path,names,ids,counts[indices])
        passed=all(metrics[k]>=v for k,v in config['minimum'].items())
        curve.append({'features':size,'bytes':path.stat().st_size,'size_ratio':path.stat().st_size/source_size,'metrics':metrics,'pass':passed})
        if selected is None and passed:
            selected={'features':size,'metrics':metrics,'bytes':path.stat().st_size,'indices':indices}
    report={'dataset_id':model.id,'status':'diagnostic_only','source_features':len(genes),'source_samples':len(ids),'source_bytes':source_size,'source_sha256':digest(source/'pasilla_gene_counts.tsv'),'sample_order':ids,'baseline':config['baseline'],'contrast':full['contrast'],'candidate_curve':curve,'smallest_passing_features':None if selected is None else selected['features'],'rights':'AMBER','release_eligible':False,'versions':{'python':platform.python_version(),'numpy':np.__version__},'config_sha256':digest(ROOT/'workflows/pasilla-source.json'),'workflow_sha256':digest(Path(__file__))}
    write_json(staging/'diagnostic-report.json',report)
    if selected:
        idx=selected['indices'];pocket=diagnostic(counts[idx],model.samples)
        write_json(staging/'plot-data.json',{'sample_ids':ids,'conditions':[s.condition for s in model.samples],'full_pca':full['pca'].tolist(),'pocket_pca':pocket['pca'].tolist(),'full_explained':full['explained'].tolist(),'pocket_explained':pocket['explained'].tolist(),'full_effects':full['effects'][idx].tolist(),'pocket_effects':pocket['effects'].tolist(),'genes':[genes[i] for i in idx]})
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--staging',type=Path,required=True,help='Local working directory, outside the published datasets tree.')
    args=p.parse_args()
    if args.staging.resolve().is_relative_to((ROOT/'datasets').resolve()):p.error('Use a staging directory outside datasets; rights review is pending.')
    print(json.dumps(run(args.staging.resolve()),indent=2))
