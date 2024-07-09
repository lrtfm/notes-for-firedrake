from petsc4py import PETSc

class YAOptionManager(object):
    commandline_options = frozenset(PETSc.Options().getAll())
    
    def __init__(self, parameters, options_prefix=None, cmd_high_priority=True):
        """Yet Another OptionsManager
        
        Args:
            parameters: a dict of parameters
            options_prefix: options prefix, default: None, i.e ''
            cmd_high_priority: commandline args with highest priority if true, else with lowest priority.
        """

        if options_prefix == None:
            options_prefix = '' 

        if len(options_prefix) and not options_prefix.endswith("_"):
            options_prefix += "_"
        
        self.options_prefix = options_prefix

        if parameters == None:
            parameters = {}
        else:
            parameters = flatten_parameters(parameters)
        if cmd_high_priority:
            parameters = {k: v for k, v in parameters.items() 
                          if self.options_prefix + k not in self.commandline_options}

        self.parameters = parameters 

    def __enter__(self):
        self.opts = PETSc.Options()
        self.to_restore = {}
        self.to_delete = {}
        for k, v in self.parameters.items():
            key = self.options_prefix + k
            if key in self.opts: 
                if v != self.opts[key]:
                    self.to_restore[key] = self.opts[key]
            else:
                self.to_delete[key] = v
            self.opts[key] = v

    def __exit__(self ,type, value, traceback): 
        for k in self.to_delete.keys():
            self.opts.delValue(k)
        for k, v in self.to_restore.items():
            self.opts[k] = v

