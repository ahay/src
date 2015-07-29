class Parameter():
    '''
    A Parameter contains all of the information about a Parameter
    that a Program requires (or uses) on the command line.
    
    Ideally, this would be the rsfpar class from rsf.doc, but again
    this class is shorter, and cleaner than rsfpar is.
    '''
    
    def __init__(self,name,type,desc,default,range):
        
        self.name = name
        self.desc = desc
        self.type = type
        self.range = range
        self.default = default      
        self.value = None
        
    def set(self,value):
        '''
        Check and set the value for this parameter.
        value should be a string when it comes into this function.
        '''
        
        if not value or value == '':
            self.value = value
            return
            
        if '#' in self.name:
            for tval in value.split(' '):
                if '=' in tval:
                    tval = tval.split('=')[1]
                else:
                    raise Exception("%s must specify which axis to operate on (1-9).  Example: %s=%s" % (self.name, self.name.replace('#','1'),'1'))
                if not self.check(tval):
                    raise Exception("Bad value for this parameter")
                
            self.value = value 
        elif self.check(value):
            self.value = value
        else:
            raise Exception("%s is not a valid value for this parameter." % value)
        
    def check(self,value):
        '''
        Determine if this value is a valid option for this parameter.
        '''
        valid = True
        if self.type == 'bool':
            if value == '1' or value == '0' or value == 'y' or value == 'n':
                valid = True
            else:
                valid = False
        elif self.type == 'int':
            try:
                val = int(value)
            except:
                valid = False
        elif self.type == 'float':
            try:
                val = float(value)
            except:
                valid = False
        elif self.type == 'ints':
            try:
                for item in value.split(','):
                    int(item)
            except:
                valid = False
        elif self.type == 'floats':
            try:
                for item in value.split(','):
                    float(item)
            except:
                valid = False
        elif self.type == 'bools':
            try:
                for item in value.split(','):
                    if item == '1' or item == '0' or item == 'y' or item == 'n':
                        pass
                    else:
                        raise Exception('bad value')
            except:
                valid = False
        else: # String?
            valid = True
        
        return valid

    def __str__(self):
        return '(%s) %s=%s %s - %s' % (self.type,self.name,
                                       self.default,self.range,self.desc)
                                       
    def copy(self):
        ''' 
        Returns a deep copy of this Parameter.
        '''
        name = self.name
        desc = self.desc
        type = self.type
        range = self.range
        default = self.default
        return Parameter(name,type,desc,default,range)
        


