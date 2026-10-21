from pathlib import Path
import os
import nbformat
from IPython.core.interactiveshell import InteractiveShell
from IPython.utils.capture import capture_output
import matplotlib
matplotlib.use('module://matplotlib_inline.backend_inline')
p=Path(__file__).parent/'Quantum_Splay_Certification.ipynb'
os.chdir(p.parent)
n=nbformat.read(p,as_version=4)
shell=InteractiveShell.instance()
from matplotlib_inline.backend_inline import configure_inline_support
configure_inline_support(shell,'module://matplotlib_inline.backend_inline')
count=0
for index,cell in enumerate(n.cells):
 if cell.cell_type!='code':continue
 count+=1;print('Executing cell',index,flush=True)
 cell.outputs=[];cell.execution_count=count
 with capture_output() as captured:
  result=shell.run_cell(cell.source,store_history=False)
 if captured.stdout: cell.outputs.append(nbformat.v4.new_output('stream',name='stdout',text=captured.stdout))
 if captured.stderr: cell.outputs.append(nbformat.v4.new_output('stream',name='stderr',text=captured.stderr))
 for out in captured.outputs:
  cell.outputs.append(nbformat.v4.new_output('display_data',data=out.data,metadata=out.metadata))
 nbformat.write(n,p)
 if result.error_before_exec or result.error_in_exec:
  raise RuntimeError(str(result.error_before_exec or result.error_in_exec))
print('Completed',count,'code cells',flush=True)
