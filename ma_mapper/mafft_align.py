import subprocess

def mafft_wrapper(input_filepath, nthread = None, nthreadtb = None, output_filepath = None, mafft_arg = '--threadit 0 --treeout --reorder '):
    mafft_command = 'mafft '
    if nthread is not None:
        nthread_arg = '--thread '+ str(nthread) + ' ' 
        mafft_command += nthread_arg
    if nthreadtb is not None:
        nthreadb_arg = '--threadb '+ str(nthread) + ' ' 
        mafft_command += nthreadb_arg
    if mafft_arg is not None:
        mafft_command += mafft_arg
    mafft_command += input_filepath
    if output_filepath is None:
        output_filepath = input_filepath+'.aligned' 
    try:
        output = subprocess.check_output(mafft_command, shell = True)
        with open(output_filepath, 'wb') as file:
            file.write(output)
    except subprocess.CalledProcessError as e:
        print(e.cmd)
        print(e.output)
        print(e.returncode)
        print("MAFFT error")
        raise