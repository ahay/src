from rsf.tex import *

if 'animate.pdf':
    os.environ['PSTEXPENOPTS']='serifs=n fat=2.0 fatmult=2.0 txscale=1.5 color=y force=y'

Paper('animate',lclass='beamer',options='mathserif',
      use='listings,color,amsmath,amsfonts,helvet,animate',
      include=r'''
      \mode<presentation>
      {
      \beamertemplatenavigationsymbolsempty
      \usetheme{default}
      \useoutertheme{default}
      % Background gradient and colors
      \setbeamercolor{structure}{bg=black, fg=white}
      \setbeamertemplate{frametitle}[default][center]
      \setbeamercolor{normal text}{bg=black, fg=white}
      \setbeamercolor{section in toc}{use=structure,parent=structure,fg=structure.fg!50}
      % Elements style
      \setbeamertemplate{itemize items}[circle]
      \setbeamertemplate{enumerate items}[default]
      % Font sizes
      \setbeamerfont{frametitle}{size=\LARGE}
      \setbeamerfont{title}{size=\huge}
      % Make structure elements bold
      \setbeamerfont{normal text}{series=\bfseries}
      \setbeamerfont{alerted text}{series=\bfseries}
      \setbeamerfont{example text}{series=\bfseries}
      \setbeamerfont{block body}{series=\bfseries}
      \setbeamerfont{caption}{series=\bfseries}
      \setbeamerfont{date}{series=\bfseries}
      \setbeamerfont{institute}{series=\bfseries}
      \setbeamerfont{author}{series=\bfseries}
      \setbeamerfont{itemize/enumerate body}{series=\bfseries} 
      \usefonttheme{structurebold}
      }
      ''',
      color='ALL',hires='ALL')

End()

