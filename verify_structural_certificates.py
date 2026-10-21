from structural_model import *
from pathlib import Path
import json,sys
if __name__=='__main__':
    directory=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).parent/'splay_outputs'
    paths=sorted(directory.glob('structural_*_n*.json'))
    if not paths:raise SystemExit('No structural outcomes found')
    graphs={}
    for path in paths:
        payload=json.loads(path.read_text());n=payload['n']
        if n not in graphs:graphs[n]=build_graph(n)
        print(verify_structural_payload(payload,graphs[n]),flush=True)
