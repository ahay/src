from Flow import Flow
from Program import plotPrograms 
from Program import plotCommands
from Program import Program
from Parameter import Parameter


def checkNameType(prog,ftype):
    names = ['sf'+name for name in plotPrograms.split(' ')]
    if not prog in names:
        return True
        
    if (ftype == Flow.PLOT or ftype == Flow.RESULT):
        return True
    else:
        return False
        

def getPrograms():
    ''' Get a dictionary of Programs and a list of Program names from the
    RSF self-documentation system.
    '''
    import rsf.doc
    import rsf.prog
    
    programs = {}
    for progname in rsf.doc.progs.keys():
        
        prog = rsf.doc.progs[progname]
        
        name     = prog.name
        synopsis = prog.snps
        comments = prog.cmts
        desc     = prog.desc
        file     = prog.file
        also     = prog.also
        wiki     = prog.wiki
        vers     = prog.vers
        uses     = prog.uses.copy()
        program = Program(name,synopsis,comments,desc,uses,also,wiki,file,vers)
        
        for parname,par in prog.pars.items():
            
            type = par.type.replace('\t','').replace('\n','').replace(' ','')
            desc = par.desc.replace('\t','').replace('\n','')
            default = par.default.replace('\t','').replace('\n','').replace(' ','')
            if default == '=': default = ''
            prange = par.range.replace('\t','').replace('\n','')
           
            parameter = Parameter(parname,type,desc,default,prange)
            program.add(parameter)
            
        custom = Parameter('custom','string','additional parameters specified by the user without key=val formatting','','')
        input = Parameter('input','string','input rsf file names','','')
        output = Parameter('output','string','output rsf file names','','')
        program.add(custom)
        program.add(input)
        program.add(output)
        programs[progname] = program
    
    
    for plotcommand in plotCommands.split(' '):
        name = plotcommand
        synopsis = 'Vppen command'
        comments = ''
        desc = ''
        file = ''
        also = ''
        wiki = ''
        vers = ''
        uses = {}
        program = Program(name,synopsis,comments,desc,uses,also,wiki,file,vers)
        custom = Parameter('custom','string','additional parameters specified by the user without key=val formatting','','')
        input = Parameter('input','string','input rsf file names','','')
        output = Parameter('output','string','output rsf file names','','')
        program.add(custom)
        program.add(input)
        program.add(output)
        programs[plotcommand] = program
        
    programNames = programs.keys()
    programNames.sort()

    return programs, programNames
