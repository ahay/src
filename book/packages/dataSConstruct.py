#########################################
#                                       #
# ##################################### #
# #  Data Processing Script for       # #
# #         rsf/book/data             # #
# ##################################### #
#                                       #
#########################################

from rsfproj import *
class process:
    def __init__(self,input):   
  ####################################
  #          Input Data              #
  ####################################
        for file in input['segyFiles']:
            Fetch(file,input['source'])
        for file in input['nativeFiles']:
            Fetch(file,input['source'])
  ####################################
  #   Convert Data to RSF format     #
  ####################################
        RSFfiles=[]
  #------- Convert SEGY Files -------#   
        segyHeaderNumber=0
        for file in input['segyFiles']:
            fileOut = ''
            fileLength =len(file)
            for letter in file[0:fileLength]:
                if letter is not '.':
                    fileOut = fileOut + letter
                if letter == '.':
                    break
            if input['segyHeader'](segyHeaderNumber) is not '':
                Flow(fileOut,file,'segyread tape=$SOURCE | put "segyHeader[segyHeaderNumber]"',stdin=0)
                segyHeaderNumber=segyHeaderNumber+1
                RSFfiles.append(fileOut)
            else:
                Flow(fileOut,file, 'segyread tape=$SOURCE',stdin=0)
                RSFfiles.append(fileOut)
    #----- Convert Native Files -----#
        nativeHeaderNumber=0
        for file in input['nativeFiles']:
            fileOut = ''
            fileLength = len(file)
            for letter in file[0:fileLength]:
                if letter is not '.':
                    fileOut = fileOut + letter
                if letter == '.':
                    break
            if input['nativeHeader'][nativeHeaderNumber]:
                Flow([fileOut,fileOut+'header'],file,'dd form=native | put "nativeHeader[nativeHeaderNumber]"')
                nativeHeaderNumber = nativeHeaderNumber + 1
                RSFfiles.append(fileOut)
            else:
                Flow(fileOut,file,'dd form=native')
                RSFfiles.append(fileOut)
  #######################################
  #        Concatinate Data             #
  #######################################
        if input['catRules'] is not '':
            Flow('RSFfile',RSFfiles,input['catRules'],stdin=0)
  #######################################
  #          Display Data               #
  #######################################
        if input['catRules'] is not '':
            Result('RSFfile','window j3=20 | grey gainpanel=a title="input["plotTitle"]"')
        else:
            titleNumber = 0;
            for file in RSFfiles:
                Result(file,'window j3=20 | grey gainpanel=a title=%s' % input['plotTitles'][titleNumber])
                titleNumber = titleNumber + 1
