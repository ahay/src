# coding: utf-8
#
# graph_scons.py (Python)
# 
# Purpose: Definition of ploting functions.
# 
# Site: http://www.dirackslounge.online
# 
# Version 1.0
# 
# Programer: Rodolfo Dirack 22/08/2019
# 
# Email: rodolfo_profissional@hotmail.com
# 
# License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

def wiggle(title):
    return '''
    wiggle wheretitle=top yreverse=y transp=y label1=Time unit1=s label2=CMP unit2=km pclip=99.5 title="%s" 
    ''' % title

def grey(title):
    return '''
    grey label1=Time unit1=s label2=CMP unit2=km pclip=99.5 title="%s" 
    ''' % title

def grey3(title):
	return '''
	byte |
	transp plane=23 |
	grey3 flat=n frame1=500 frame3=80 frame2=200
	label1=Time unit1=s 
	label3=Offset unit3=km 
	label2=CMP unit2=km
	title="%s" point1=0.8 point2=0.8 
	''' % title

