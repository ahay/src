plotPrograms = '''box contour contour3 dots graph graph3 grey grey3 plas pldb plotrays thplot vplotdiff wiggle'''

plotCommands = 'Overlay SideBySideIso SideBySideAniso OverUnderAniso OverUnderIso TwoRows TwoColumns Movie'

class Program():
    '''
    A Program is just a wrapper for an RSF Program, that contains
    all of the information from the RSF self-documentation construct.
    
    Ideally, this would be replaced with the rsfprog class from rsf.doc
    but that class is not designed to be used with this, and has a lot
    of extra methods that are not necessary.
    '''
    
    def __init__(self,name,synopsis,comments,desc,uses,also,wiki,file,vers):
        self.name = name
        self.synopsis = synopsis
        self.comments = comments
        self.desc     = desc
        self.pars = {}
        self.also = also
        self.wiki = wiki
        self.uses = uses
        self.file = file
        self.vers = vers
          
    def add(self,par):
        '''
        Add a Parameter instance to this program.
        '''
        self.pars[par.name] = par
    
    def selfdoc(self):
        '''
        Convert this Program to a self-documentation string.
        '''
        desc = '\n[NAME]:\n\t%s - %s \n\n' % (self.name,self.desc)
        snps = '[SYNOPSIS]:\n\t%s \n\n' % self.synopsis
        cmts = '[COMMENTS]:\n\t%s\n\n' % self.comments
        pars = '[PARAMETERS]:\n\n'
        names = self.pars.keys()
        names.sort()
        for name in names:
            par = self.pars[name]
            pars += '\t'+str(par)+'\n'
        uses = '\n[USED IN]:\n'
        books = self.uses.keys()
        books.sort()
        for book in self.uses.keys():
            chapters = self.uses[book].keys()
            chapters.sort()
            for chapter in chapters:
                for project in self.uses[book][chapter]:
                    uses += '\t%s/%s/%s\n' % (book,chapter,project)
        source = '\n[SOURCE]:\n\t%s\n' % self.file
        wiki   = '\n[WIKI]:\n\t%s\n' % self.wiki
        if not self.wiki or self.wiki == '': wiki = ''
        vers   = '\n[VERSION]:\n\t%s\n' % self.vers
        if not self.vers or self.vers == '': vers = ''
        also   = '\n[SEE ALSO]:\n\t%s\n' % self.also
        if not self.also or self.also == '': also = ''
        return desc+snps+cmts+pars+uses+source+wiki+vers+also
        
    def __str__(self):
        return self.selfdoc()
