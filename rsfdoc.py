import pydoc

def bold(text):
    """Format a string in bold by overstriking."""
    return ''.join(map(lambda ch: ch + "\b" + ch, text))

def show(head,body):
    text = "\n".join(map(lambda line: "\t" + line, body.split("\n")))
    pydoc.pager(bold(head.upper()) + "\n" + text)

class par:
    def __init__(self,name,type,default=None,desc=None):
        self.name = name
        self.type = type
        self.default = self.default
        self.desc = self.desc
    
class rsfprog:
    def __init__(self,name,file,desc=None):
        self.name = name
        self.file = file
        self.desc = desc
        self.pars = {}
    def add_par (par):
        self.pars[par.name] = par
    def doc(self):
        show('name',self.name)
        if self.desc:
            show('description',self.desc)

sfin = rsfprog('sfin','sfin.c')
sfin.doc()
