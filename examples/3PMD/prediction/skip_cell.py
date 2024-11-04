def run_if(line, cell=None):
    '''Runs the cell with the skip magic.'''
    if eval(line) == False:
        print("[run_if] Skipping cell")
        return 
    else:
        get_ipython().ex(cell)

def skip_if(line, cell=None):
    '''Skips execution of the current line/cell if line evaluates to True.'''
    if eval(line) == True:
        print("[skip_if] Skipping cell")
        return
    else:
        get_ipython().ex(cell)

def load_ipython_extension(shell):
    '''Registers the skip magic when the extension loads.'''
    shell.register_magic_function(skip_if, 'line_cell')
    shell.register_magic_function(run_if, 'line_cell')

def unload_ipython_extension(shell):
    '''Unregisters the skip magic when the extension unloads.'''
    del shell.magics_manager.magics['cell']['skip_if']
    del shell.magics_manager.magics['cell']['run_if']