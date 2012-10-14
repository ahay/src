%%
%% This is file `geophysics.cls',
%% generated with the docstrip utility.
%%
%% The original source files were:
%%
%% geophysics.dtx  (with options: `class')
%% 
%% This is an automatically generated file. Do not edit it directly!
%% 
%%    Copyright (C) 2004 by Sergey Fomel
%% 
%%     This work may be distributed and/or modified under the
%%     conditions of the LaTeX Project Public License, either version 1.3
%%     of this license or (at your option) any later version.
%%     The latest version of this license is in
%%       http://www.latex-project.org/lppl.txt
%%     and version 1.3 or later is part of all distributions of LaTeX
%%     version 2003/12/01 or later.
%% 
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{geophysics}[2007/02/14 v3.0 Geophysics paper]
%% \CharacterTable
%%  {Upper-case    \A\B\C\D\E\F\G\H\I\J\K\L\M\N\O\P\Q\R\S\T\U\V\W\X\Y\Z
%%   Lower-case    \a\b\c\d\e\f\g\h\i\j\k\l\m\n\o\p\q\r\s\t\u\v\w\x\y\z
%%   Digits        \0\1\2\3\4\5\6\7\8\9
%%   Exclamation   \!     Double quote  \"     Hash (number) \#
%%   Dollar        \$     Percent       \%     Ampersand     \&
%%   Acute accent  \'     Left paren    \(     Right paren   \)
%%   Asterisk      \*     Plus          \+     Comma         \,
%%   Minus         \-     Point         \.     Solidus       \/
%%   Colon         \:     Semicolon     \;     Less than     \<
%%   Equals        \=     Greater than  \>     Question mark \?
%%   Commercial at \@     Left bracket  \[     Backslash     \\
%%   Right bracket \]     Circumflex    \^     Underscore    \_
%%   Grave accent  \`     Left brace    \{     Vertical bar  \|
%%   Right brace   \}     Tilde         \~}
\RequirePackage{ifthen}
\newboolean{@manu}
\setboolean{@manu}{false}
\DeclareOption{manuscript}{\setboolean{@manu}{true}}
\DeclareOption{paper}{\setboolean{@manu}{false}}
\newboolean{@twoc}
\setboolean{@twoc}{false}
\DeclareOption{twocolumn}{\setboolean{@twoc}{true}}
\DeclareOption{onecolumn}{\setboolean{@twoc}{false}}
\newboolean{@shrt}
\setboolean{@shrt}{false}
\DeclareOption{short}{\setboolean{@shrt}{true}}
\DeclareOption{long}{\setboolean{@shrt}{false}}
\newcommand{\shortpaper}{\setboolean{@shrt}{true}}
\newcommand{\longpaper}{\setboolean{@shrt}{false}}
\newboolean{@eflt}
\setboolean{@eflt}{false}
\DeclareOption{endfloat}{\setboolean{@eflt}{true}}
\DeclareOption{noendfloat}{\setboolean{@eflt}{false}}
\newboolean{@revd}
\setboolean{@revd}{false}
\DeclareOption{revised}{\setboolean{@revd}{true}}
\DeclareOption{new}{\setboolean{@revd}{false}}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ExecuteOptions{paper,onecolumn,long,noendfloat,new}
\ProcessOptions*
\ifthenelse{\boolean{@twoc}}{%
  \LoadClass[9pt]{article}}{%
  \LoadClass[11pt]{article}}
\RequirePackage{seg}
\ifthenelse{\boolean{@eflt}}{
  \ifthenelse{\boolean{@manu}}{\RequirePackage{endfloat}}{}
}{
\ifthenelse{\boolean{@manu}}{
  \renewcommand{\inputdir}[1]{\AtEndDocument{\renewcommand{\@path}{#1}}}}{}
\providecommand{\figdir}{.}
\providecommand{\setfigdir}[1]{\renewcommand{\figdir}{#1}}
\providecommand{\fig@file}[1]{\@path/\figdir/#1}
\newcounter{@plots}
\setcounter{@plots}{0}
\ifthenelse{\boolean{@manu}}
{\newcommand{\fullplot}[4][X]{
    \stepcounter{@plots}
    \AtEndDocument{
      \begin{figure}[p]
        \centering
        \includegraphics[#3]{\fig@file{#2}}
        \caption [#4]{#4 \\
  \textbf{\seg@lhead} --
  \ifthenelse{\equal{\ms@number}{Manuscript}}{}{\ms@number}}
\label{fig:\@path-#2}
        \label{fig:#2}
      \end{figure}\clearpage
    }}
  \def\plot{\@ifstar{\fullplot}{\fullplot}}
  \setlength\abovecaptionskip{50\p@}
}
{\newcommand{\fullplot}[4][htbp]{%
    \begin{figure*}[#1]
      \centering
      \includegraphics[#3]{\fig@file{#2}}
      \caption{#4}
      \label{fig:\@path-#2}
      \label{fig:#2}
    \end{figure*}
  }}
\ifthenelse{\boolean{@manu}}
{\newcommand{\sideplot}[4][X]{\plot[#1]{#2}{#3}{#4}}}
{\newcommand{\sideplot}[4][htbp]{%
    \ifthenelse{\boolean{@twoc}}{%
      \begin{figure}[#1]
        \centering
        \includegraphics[#3]{\fig@file{#2}}
        \caption{#4}
\label{fig:\@path-#2}
        \label{fig:#2}
      \end{figure}
    }{%
      \begin{figure*}[#1]
        \makebox{%
          \begin{minipage}{0.4\textwidth}
            \caption{#4}
    \label{fig:\@path-#2}
            \label{fig:#2}
          \end{minipage}
          \hfill
          \begin{minipage}{0.6\textwidth}
            \centering
            \includegraphics[#3]{\fig@file{#2}}
          \end{minipage}}
      \end{figure*}
    }
  }
  \ifthenelse{\boolean{@twoc}}{%
    \def\plot{\@ifstar{\fullplot}{\sideplot}}}{%
    \def\plot{\@ifstar{\fullplot}{\fullplot}}}
}
\RequirePackage{subfigure}
\def\next@item#1,#2?{#1}
\def\rest@item#1,#2?{#2}
\newcounter{sub@fig}
\ifthenelse{\boolean{@manu}}
{\newcommand{\seg@multiplot}[5][X]{
    \stepcounter{@plots}
    \AtEndDocument{
      \begin{figure}[p]
        \centering
        \setcounter{sub@fig}{0}
        \edef\list@i{#3}
        \loop
        \edef\item@i{\expandafter\next@item\list@i,\empty?}
        \edef\list@i{\expandafter\rest@item\list@i,\empty?}
        \ifx\item@i\empty
        \else
        \stepcounter{sub@fig}
        \subfigure[]{\includegraphics[#4]{\fig@file{\item@i}}%
  \label{fig:\@path-\item@i}
  \label{fig:\item@i}}
        \ifthenelse{\value{sub@fig} = #2}{\\ \setcounter{sub@fig}{0}}{}
        \repeat
        \caption[#5]{#5 \\
  \textbf{\seg@lhead} --
  \ifthenelse{\equal{\ms@number}{Manuscript}}{}{\ms@number}}
\label{fig:\@path-#3}
        \label{fig:#3}
      \end{figure}\clearpage
}}
\def\multiplot{\@ifstar{\seg@multiplot}{\seg@multiplot}}
}{%
  \newcommand{\seg@fullmultiplot}[5][htbp]{
    \begin{figure*}[#1]
      \centering
      \setcounter{sub@fig}{0}
      \edef\list@i{#3}
      \loop
      \edef\item@i{\expandafter\next@item\list@i,\empty?}
      \edef\list@i{\expandafter\rest@item\list@i,\empty?}
      \ifx\item@i\empty
      \else
      \stepcounter{sub@fig}
      \subfigure[]{\includegraphics[#4]{\fig@file{\item@i}}%
\label{fig:\@path-\item@i}
        \label{fig:\item@i}}
      \ifthenelse{\value{sub@fig} = #2}{\\ \setcounter{sub@fig}{0}}{}
      \repeat
      \caption{#5}
      \label{fig:\@path-#3}
      \label{fig:#3}
    \end{figure*}
  }
  \newcommand{\seg@sidemultiplot}[5][htbp]{
    \begin{figure}[#1]
      \centering
      \setcounter{sub@fig}{0}
      \edef\list@i{#3}
      \loop
      \edef\item@i{\expandafter\next@item\list@i,\empty?}
      \edef\list@i{\expandafter\rest@item\list@i,\empty?}
      \ifx\item@i\empty
      \else
      \stepcounter{sub@fig}
      \subfigure[]{\includegraphics[#4]{\fig@file{\item@i}}%
\label{fig:\@path-\item@i}
        \label{fig:\item@i}}
      \ifthenelse{\value{sub@fig} = #2}{\\ \setcounter{sub@fig}{0}}{}
      \repeat
      \caption{#5}
      \label{fig:\@path-#3}
      \label{fig:#3}
    \end{figure}
  }
  \ifthenelse{\boolean{@twoc}}{%
    \def\multiplot{\@ifstar{\seg@fullmultiplot}{\seg@sidemultiplot}}}{%
    \def\multiplot{\@ifstar{\seg@fullmultiplot}{\seg@fullmultiplot}}}
}
\newcounter{@tabls}
\setcounter{@tabls}{0}
\ifthenelse{\boolean{@manu}}{
  \providecommand{\tabl}[4][X]{
    \stepcounter{@tabls}
    \AtEndDocument{
      \begin{table}[p]
        #4
        \caption{#3}
        \label{tbl:#2}
      \end{table}\clearpage
    }
  }
}{
  \providecommand{\seg@sidetabl}[4][htbp]{
    \begin{table}[#1]
      #4
      \caption{#3}
      \label{tbl:#2}
    \end{table}
  }
  \providecommand{\seg@fulltabl}[4][htbp]{
    \begin{table*}[#1]
      #4
      \caption{#3}
      \label{tbl:#2}
    \end{table*}
  }
  \def\tabl{\@ifstar{\seg@fulltabl}{\seg@sidetabl}}
}
\renewcommand*{\l@figure}[2]{%
  \setlength{\@tempdima}{2.3em}%
  \noindent\hspace*{1.5em}#1\hfil\newline}
\let\l@table\l@figure
\renewcommand\listoffigures{%
  \section{\listfigurename
    \@mkboth{\MakeUppercase\listfigurename}%
    {\MakeUppercase\listfigurename}}%
  \@starttoc{lof}%
}
\renewcommand\listoftables{%
  \section{\listtablename
    \@mkboth{\MakeUppercase\listtablename}%
    {\MakeUppercase\listtablename}}%
  \@starttoc{lot}%
}
\AtEndDocument{
  \ifthenelse{\value{@tabls} > 0}{\newpage\listoftables}{}
  \ifthenelse{\value{@plots} > 0}{\newpage\listoffigures}{}
}}
\ifthenelse{\boolean{@revd}}{
  \RequirePackage{color}%
  \RequirePackage{ulem}%
  \normalem
}{}
\DeclareRobustCommand{\new}[1]{%
  \ifthenelse{\boolean{@revd}}{\textit{\textcolor{red}{#1}}}{#1}}
\DeclareRobustCommand{\old}[1]{\ifthenelse{\boolean{@revd}}{\textcolor{blue}{\sout{#1}}}{}}
\RequirePackage{natbib}
\bibpunct[,]{(}{)}{;}{a}{,}{,}
\DeclareRobustCommand{\refname}{REFERENCES}
\ifthenelse{\boolean{@manu}}
{\renewcommand\bibsection{\newpage\section{\refname}}}
{\renewcommand\bibsection{\section{\refname}}}
\newcommand{\shortcite}[1]{\citeyearpar{#1}}
\ifthenelse{\boolean{@manu}}{
  \renewcommand{\baselinestretch}{2.0}
}{
  \renewcommand{\baselinestretch}{1.0}
}
\ifthenelse{\boolean{@twoc}}{%
\setlength\parindent{1em}
\setlength\headheight{12.5\p@}
\setlength\headsep   {11.9\p@}
\setlength\topskip   {3\p@}
\setlength\footskip{19pt}
\setlength\maxdepth{4\p@}
\setlength\@maxdepth\maxdepth
\setlength\textwidth{42pc}
\addtolength\textwidth{2em}
\newdimen\typeheight \typeheight55.5pc
\advance\typeheight-2.215pt
\setlength\textheight{\typeheight}
\addtolength\textheight{\headheight}
\addtolength\textheight{\headsep}
\setlength\columnsep{24pt}
\setlength\columnseprule{0\p@}
\setlength\oddsidemargin{-18pt}
\setlength\evensidemargin{-18pt}
\addtolength\oddsidemargin{-1em}
\addtolength\evensidemargin{-1em}
\setlength\marginparwidth {.75in}
\setlength\marginparsep {\z@}
\setlength\marginparpush{5\p@}
\setlength\topmargin{-29.5pt}
\addtolength{\topmargin}{-\headheight}
\addtolength{\topmargin}{-\headsep}
\setlength\footnotesep{5.65\p@}%{4pt}%
\setlength{\skip\footins}{12\p@ \@plus 2\p@}%
\setlength\floatsep    {16\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\textfloatsep{20\p@ \@plus 2\p@ \@minus 4\p@}
\setlength\intextsep   {12\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\dblfloatsep    {12\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\dbltextfloatsep{20\p@ \@plus 2\p@ \@minus 4\p@}
\setlength\@fptop{0\p@ \@plus 1fil}
\setlength\@fpsep{8\p@ \@plus 2fil}
\setlength\@fpbot{0\p@ \@plus 1fil}
\setlength\@dblfptop{0\p@ \@plus 1fil}
\setlength\@dblfpsep{8\p@ \@plus 2fil}
\setlength\@dblfpbot{0\p@ \@plus 1fil}
\setlength\partopsep{2\p@ \@plus 1\p@ \@minus 1\p@}
\setlength\lineskip{\z@}%
\setlength\normallineskip{\z@}%
\setlength\parskip{0\p@}
\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301
\setlength\leftmargini  {2.4em}
\setlength\leftmarginii  {2.2em}
\setlength\leftmarginiii {1.87em}
\setlength\leftmarginiv  {1.7em}
\setlength\leftmarginv  {.5em}
\setlength\leftmarginvi {.5em}
\setlength\leftmargin    {\leftmargini}
\setlength  \labelsep  {.5em}
\setlength  \labelwidth{\leftmargini}
\addtolength\labelwidth{-\labelsep}
\@beginparpenalty -\@lowpenalty
\@endparpenalty   -\@lowpenalty
\@itempenalty     -\@lowpenalty
}{%
\setlength{\textheight}{9in}
\setlength{\textwidth}{6in}
\setlength{\parskip}{1.5ex plus0.1ex minus0.1ex} % should be a *rubber* length
\setlength{\oddsidemargin}{.26in}
\setlength{\evensidemargin}{.26in}
\setlength{\topmargin}{-.51in}
\setlength{\textfloatsep}{30pt plus 3pt minus 6pt}
\setlength{\headsep}{0.5in}}
\let\seg@large\large
\let\seg@Large\Large
\renewcommand{\Large}{\protect\seg@large}
\let\seg@section\section
\let\seg@subsection\subsection
\newcommand{\segsection}{%
\@startsection {section}{1}{\z@}%
{-3.5ex \@plus -1ex \@minus -.2ex}%
{1ex \@plus .2ex}%
{\centering%
\normalfont\ifthenelse{\boolean{@twoc}}{}{\Large}\bfseries\MakeUppercase}}
\newcommand{\segsubsubsection}{%
\@startsection{subsubsection}{3}{\z@}%
{-3.25ex\@plus -1ex \@minus -.2ex}%
{1.5ex \@plus .2ex}%
{\normalfont\normalsize\itshape}}
\newcommand{\seg@subsubsubsection}[2][]{%
\underline{#2}.--}
\def\section{\@ifstar{\segsection*}{\segsection*}}
\def\subsection{\@ifstar{\seg@subsection*}{\seg@subsection*}}
\def\subsubsection{\@ifstar{\segsubsubsection*}{\segsubsubsection*}}
\def\subsubsubsection{\@ifstar{\seg@subsubsubsection}{\seg@subsubsubsection}}
\newcommand{\ms@number}{Manuscript}
\newcommand{\ms}[1]{\renewcommand{\ms@number}{\textbf{#1}}}
\newcommand{\geo@published}{Unknown}
\newcommand{\published}[1]{\renewcommand{\geo@published}{#1}}
\newcommand{\seg@title}{}
\newenvironment{geo@title}{%
  \begin{minipage}{\textwidth}%
    \bfseries\seg@Large\ifthenelse{\boolean{@twoc}}{\flushleft}{\centering}}
       {\end{minipage}\vspace{.2in}}%
\renewcommand{\title}[1]{%
  \renewcommand{\seg@title}{%
    \noindent\begin{geo@title}\setlength{\parindent}{0em}\par%
   {\smallskip #1\ifthenelse{\equal{\geo@published}{Unknown}}{}{\footnote{Published in \geo@published}}}\end{geo@title}}}
\providecommand{\seg@address}{}
\providecommand{\address}[1]{\renewcommand{\seg@address}{#1}}
\providecommand{\seg@author}{}
\renewcommand{\author}[1]{%
  \renewcommand{\seg@author}{%
    \ifthenelse{\boolean{@twoc}}{%
      \seg@title \\
      \begin{flushleft}
\Large #1
      \end{flushleft}\vspace{0.2in}
    }{%
      \seg@title
      \begin{center}
        \ifthenelse{\boolean{@manu}}{%
          \textbf{#1} \\
          \emph{\seg@address}\ \\
          (\today) \\
  \ifthenelse{\equal{\ms@number}{Manuscript}}{}{\ms@number \\}
  Running head: \textbf{\seg@rhead} \\
        }{\textit{#1}\ifthenelse{\equal{\seg@email}{.}}{}{%
            \footnote{\textbf{e-mail: }\seg@email}} \\ \
        }
      \end{center}
    }}}
\renewcommand{\maketitle}{\ifthenelse{\boolean{@shrt}}{%
    \ifthenelse{\boolean{@twoc}}{\twocolumn[\seg@author]}{\seg@author}}{}}
\ifthenelse{\boolean{@manu}}{
  \renewenvironment{abstract}{\seg@author
    \section{ABSTRACT}}{\newpage}
}
{
  \newsavebox{\seg@abstract}
  \setlength{\fboxsep}{10pt}
  \renewenvironment{abstract}{
    \ifthenelse{\boolean{@twoc}}{%
      \twocolumn[\vspace{1.8in}\seg@author]\footnotetext{\protect\seg@address}%
    }{\vspace{0.3in}\seg@author}
    \begin{lrbox}{\seg@abstract}\begin{minipage}{0.95\columnwidth}%
         \centerline{\ifthenelse{\boolean{@twoc}}{}{\Large}\textbf{ABSTRACT}}\medskip}{%
    \medskip\end{minipage}\end{lrbox}%
    \ifthenelse{\boolean{@shrt}}{}{\noindent\fbox{\usebox{\seg@abstract}}}}
}
\ifthenelse{\boolean{@manu}}
{\pagestyle{plain}\thispagestyle{plain}}
{\newcommand{\ps@seg}{%
\ifthenelse{\boolean{@twoc}}{%
\renewcommand{\@oddhead}{\hfill\textbf{\seg@rhead}\hfill\textbf{\thepage}}
\renewcommand{\@evenhead}{\textbf{\thepage}\hfill\textbf{\seg@lhead}\hfill}
\renewcommand{\@oddfoot}{}
\renewcommand{\@evenfoot}{}}{%
\renewcommand{\@oddhead}{%
\makebox[2.5in][l]{\textit{\seg@lhead}}\hfill\thepage\hfill%
\makebox[2.5in][r]{\textit{\seg@rhead}}}%
\renewcommand{\@evenhead}{\@oddhead}%
\renewcommand{\@oddfoot}{\hfill\textit{\seg@foot}\hfill}%
\renewcommand{\@evenfoot}{\@oddfoot}}}
\newcommand{\ps@segone}{%
\renewcommand{\@oddhead}{}%
\renewcommand{\@evenhead}{}%
\ifthenelse{\boolean{@twoc}}{%
  \renewcommand{\@oddfoot}{\hfill\textbf{\thepage}\hfill}
}{\renewcommand{\@oddfoot}{\hfill\textit{\seg@foot}\hfill}}%
\renewcommand{\@evenfoot}{\@oddfoot}}
\pagestyle{seg}\thispagestyle{segone}}
\newcounter{@append}
\ifthenelse{\boolean{@manu}}{%
  \AtEndDocument{\setcounter{@append}{0}}}{}
\providecommand{\append@name}{A}
\providecommand{\appendname}[1]{\renewcommand{\append@name}{#1}}
\renewcommand{\appendix}{%
  \ifthenelse{\equal{\append@name}{A}}{%
    \renewcommand{\append@name}{\Alph{@append}}}{}
  \setcounter{equation}{0}%
  \renewcommand{\theequation}{\mbox{\append@name-\arabic{equation}}}%
  \ifthenelse{\boolean{@manu}}{%
    \AtEndDocument{
      \setcounter{figure}{0}%
      \renewcommand{\thefigure}{\append@name-\arabic{figure}}%
      \renewcommand{\p@subfigure}{\thefigure}%
      \stepcounter{@append}
    }
  }{%
    \setcounter{figure}{0}%
    \renewcommand{\thefigure}{\append@name-\arabic{figure}}%
    \renewcommand{\p@subfigure}{\thefigure}%
  }
  \stepcounter{@append}
}
\providecommand{\append}[2][\append@name]{%
  \appendix\section{Appendix \append@name}\section{#2}
  \begingroup
  \def\@currentlabel{\append@name}%
  \label{#1}%
  \endgroup
}
\endinput
%%
%% End of file `geophysics.cls'.
