"""Write the biological collection decision and fail closed when evidence is missing."""
import json
from pathlib import Path
from omicsbench.registry import registry
from omicsbench.validate import validate

root=Path(__file__).resolve().parents[1]
metadata=json.loads((root/'release/metadata-input.json').read_text())
records=registry(root)
blockers=[]
validated=[]
for model,folder in records.values():
    try:result=validate(model,folder)
    except ValueError as exc:blockers.append(str(exc));continue
    if model.kind=='real' and model.status=='validated' and result['status']=='pass':validated.append(model)
pockets=[m for m in validated if m.rights.status=='GREEN' and any(f.tier=='pocket' for f in m.files)]
if not 12<=len(validated)<=20:blockers.append(f'Biological collection: {len(validated)} certified objects; target is 12 to 20.')
if len(pockets)<8:blockers.append(f'Redistributable biological pockets: {len(pockets)}; at least 8 required.')
if len({m.archetype for m in validated})<4:blockers.append('At least four certified design archetypes are required.')
for model in validated:
    if model.derivation.commit is None:blockers.append(f'{model.id}: transformation commit is missing.')
    if model.validation.baseline_version.startswith('diagnostic'):blockers.append(f'{model.id}: diagnostic-only baseline is not biological certification.')
for field in ['authors','repository_url','code_license','metadata_license','version_doi','external_user_trial','scientific_environment_lock','post_upload_verification']:
    if not metadata.get(field):blockers.append(f'Missing release evidence: {field}.')
# Positive evidence verification will be implemented with the first certified biological collection.
blockers.append('Biological collection certification is not enabled: real baseline, independent trial and deposited-file evidence adapters remain unimplemented.')
report={'scope':'planned_biological_collection','software_version':metadata['software_version'],'decision':'NO-GO' if blockers else 'GO','certified_biological_objects':len(validated),'redistributable_pockets':len(pockets),'blockers':blockers}
(root/'release/certification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(report,indent=2))
raise SystemExit(1 if blockers else 0)
