%%
%% This is file `georeport.cls',
%% generated with the docstrip utility.
%%
%% The original source files were:
%%
%% georeport.dtx  (with options: `class')
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
\ProvidesClass{georeport}[2004/11/23 v1.0 SEP-style report]
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
\newboolean{@twoc}
\setboolean{@twoc}{false}
\DeclareOption{twocolumn}{\setboolean{@twoc}{true}}
\DeclareOption{onecolumn}{\setboolean{@twoc}{false}}
\newboolean{@shrt}
\setboolean{@shrt}{false}
\newcommand{\shortpaper}{\setboolean{@shrt}{true}}
\newcommand{\longpaper}{\setboolean{@shrt}{false}}
\newboolean{@book}
\setboolean{@book}{false}
\DeclareOption{book}{\setboolean{@book}{true}}
\DeclareOption{report}{\setboolean{@book}{false}}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{book}}
\ExecuteOptions{onecolumn}
\ProcessOptions*
\ifthenelse{\boolean{@twoc}}{%
  \LoadClass[9pt]{book}}{%
  \LoadClass[11pt]{book}}
\RequirePackage{seg}
\RequirePackage{chapterbib}
\newlength{\toc@width}%
\setlength{\toc@width}{\textwidth}%
\addtolength{\toc@width}{-6em}%
\newcommand{\TOCentry}[3][.]{
\ifthenelse{\boolean{@book}}{}{%
  \noindent\hspace*{4em}{\begin{minipage}[b]{\toc@width}
      \setlength{\parindent}{-5em}%
      \ifthenelse{\equal{#1}{.}}{}{\textit {#1}, }%
      #2\dotfill\
    \end{minipage}\hfill#3}\par}}
\newcommand{\REPORT}{}
\newcommand{\GROUP}{}
\newcommand{\GEOcopyr}[2]{%
  \ \vfill
  \mainauthor{Copyright \copyright\ #1}
  \mainauthor{by #2}
  \vfill \ \newpage}
\newcommand{\geo@header}{}
\newcommand{\GEOheader}[1]{\renewcommand{\geo@header}{#1}}
\newcommand{\paper@dir}{.}
\newcommand{\paper@tag}{self}
\newcommand{\GEOpaper}[2]{%
\ifthenelse{\boolean{@book}}{}{%
  \setcounter{equation}{0}\setcounter{figure}{0}%
  \setcounter{table}{0}\setcounter{footnote}{0}%
  \renewcommand{\theequation}{\arabic{equation}}%
  \renewcommand{\thefigure}{\arabic{figure}}%
  \renewcommand{\thetable}{\arabic{table}}}%
  \cleardoublepage
  \label{\paper@tag.end}%
  \thispagestyle{segtwo}%
  \longpaper%
  \label{#2.start}
  \renewcommand{\paper@dir}{#1}
  \renewcommand{\paper@tag}{#2}
}
\RequirePackage{ifpdf}
\ifpdf
\RequirePackage[pdftex]{graphicx}
\DeclareGraphicsExtensions{.pdf,.png}
\else
\usepackage{graphicx}
\DeclareGraphicsExtensions{.eps,.ps}
\fi
\providecommand{\figdir}{.}
\providecommand{\setfigdir}[1]{\renewcommand{\figdir}{#1}}
\providecommand{\fig@file}[1]{\paper@dir/\@path/\figdir/#1}
\newcommand{\plotdir}{\figdir}
\newcommand{\plot@file}[1]{\plotdir/#1}
\newcommand{\plotbox}[2]{%
  \begin{figure*}[htbp]
    \centering
    \includegraphics[#2]{\plot@file{#1}}
\end{figure*}}
\newcommand{\fullplot}[4][htbp]{%
  \begin{figure*}[#1]
    \centering
    \setlength{\fboxsep}{0.5ex}
    \includegraphics[#3]{\fig@file{#2}}
    \caption{#4 \fbox{\paper@dir/\@path\ #2}}
    \label{fig:#2}
  \end{figure*}
}
\newcommand{\sideplot}[4][htbp]{%
  \ifthenelse{\boolean{@twoc}}{%
    \begin{figure}[#1]
      \centering
      \setlength{\fboxsep}{0.5ex}
      \includegraphics[#3]{\fig@file{#2}}
      \caption{#4 \fbox{\paper@dir/\@path\ #2}}
      \label{fig:#2}
    \end{figure}
  }{%
    \begin{figure*}[#1]
      \makebox{%
        \begin{minipage}{0.4\textwidth}
  \setlength{\fboxsep}{0.5ex}
          \caption{#4 \fbox{\paper@dir/\@path\ #2}}
          \label{fig:#2}
        \end{minipage}
        \hfill
        \begin{minipage}{0.6\textwidth}
          \flushright
          \includegraphics[#3]{\fig@file{#2}}
        \end{minipage}}
    \end{figure*}
  }
}
\ifthenelse{\boolean{@twoc}}{%
  \def\plot{\@ifstar{\fullplot}{\sideplot}}}{%
  \def\plot{\@ifstar{\fullplot}{\fullplot}}}
\RequirePackage{subfigure}
\def\next@item#1,#2?{#1}
\def\rest@item#1,#2?{#2}
\newcounter{sub@fig}
\newcommand{\seg@fullmultiplot}[5][htbp]{
  \begin{figure*}[#1]
    \centering
    \setlength{\fboxsep}{0.5ex}
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
    \caption{#5 \fbox{\paper@dir/\@path\ #3}}
    \label{fig:\@path-#3}
    \label{fig:#3}
  \end{figure*}
}
\newcommand{\seg@sidemultiplot}[5][htbp]{
  \begin{figure}[#1]
    \centering
    \setlength{\fboxsep}{0.5ex}
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
    \caption{#5 \fbox{\paper@dir/\@path\ #3}}
    \label{fig:\@path-#3}
    \label{fig:#3}
  \end{figure}
}
\ifthenelse{\boolean{@twoc}}{%
  \def\multiplot{\@ifstar{\seg@fullmultiplot}{\seg@sidemultiplot}}}{%
  \def\multiplot{\@ifstar{\seg@fullmultiplot}{\seg@fullmultiplot}}}
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
\newcommand{\new}[1]{#1}
\newcommand{\old}[1]{}
\newcommand{\ms}[1]{}
\def\first@char#1#2\end{#1}
\def\geo@input#1{\ifthenelse{\equal{\paper@dir}{.}\or\equal{\expandafter\first@char#1\end}{.}}{\@iinput{#1}}{\@iinput{\paper@dir/#1}}}
\def\input{\@ifnextchar\bgroup\geo@input\@@input}
\setlength{\textheight}{9in}
\setlength{\textwidth}{6in}
\setlength{\parskip}{1.5ex plus0.1ex minus0.1ex} % should be a *rubber* length
\setlength{\oddsidemargin}{.26in}
\setlength{\evensidemargin}{.26in}
\setlength{\topmargin}{-.51in}
\setlength{\textfloatsep}{30pt plus 3pt minus 6pt}
\setlength{\headsep}{0.5in}
\newcommand{\seg@title}{}
\newenvironment{geo@title}{%
  \begin{minipage}{\textwidth}%
    \bfseries\seg@Large\ifthenelse{\boolean{@twoc}}{\flushleft}{\centering}}
       {\end{minipage}\vspace{.2in}}%
\ifthenelse{\boolean{@book}}{}{%
\renewcommand{\title}[1]{%
  \renewcommand{\seg@title}{%
    \noindent\begin{geo@title}\setlength{\parindent}{0em}\par%
   {\smallskip #1}\end{geo@title}}}}
\newcommand{\maintitle}[1]{\ifthenelse{\boolean{@book}}{}{\title{#1}\seg@title}}
\providecommand{\seg@address}{}
\providecommand{\address}[1]{\renewcommand{\seg@address}{#1}}
\providecommand{\seg@author}{}
\newcommand{\mainauthor}[1]{
  \ifthenelse{\boolean{@twoc}}{%
    \begin{flushleft}
      \Large #1\footnote{\seg@address}
    \end{flushleft}\vspace{0.2in}
  }{%
    \begin{center}
      \textit{#1}\ifthenelse{\equal{\seg@email}{.}}{}{%
        \footnote{\textbf{e-mail: }\seg@email}} \\ \
    \end{center}
  }
}
\ifthenelse{\boolean{@book}}{}{%
\renewcommand{\author}[1]{%
  \renewcommand{\seg@author}{%
    \ifthenelse{\boolean{@twoc}}{%
      \seg@title \\
      \mainauthor{#1}
    }{%
      \seg@title
      \mainauthor{#1}
    }
  }
}}
\newcommand{\geo@ignore}[1]{}
\ifthenelse{\boolean{@book}}{
  \let\geo@toc\tableofcontents
  \renewcommand{\tableofcontents}{\geo@toc%
    \def\title{\chapter}%
    \def\author{\geo@ignore}}}{%
  \renewcommand{\maketitle}{\ifthenelse{\boolean{@shrt}}{%
      \ifthenelse{\boolean{@twoc}}{\twocolumn[\seg@author]}{\seg@author}}{}}}
\newsavebox{\seg@abstract}
\setlength{\fboxsep}{10pt}
\newenvironment{abstract}{
  \ifthenelse{\boolean{@twoc}}{%
    \twocolumn[\vspace{1.8in}\seg@author]}{\vspace{0.3in}\seg@author}
  \begin{lrbox}{\seg@abstract}\begin{minipage}{0.95\columnwidth}%
      \centerline{\ifthenelse{\boolean{@twoc}}{}{\Large}\textbf{ABSTRACT}}\medskip}{%
  \medskip\end{minipage}\end{lrbox}%
  \ifthenelse{\boolean{@shrt}}{}{\noindent\fbox{\usebox{\seg@abstract}}}}
\newenvironment{bios}{%
  \begin{list}{}{%
      \setlength{\leftmargin}{1cm}\setlength{\rightmargin}{1cm}%
      \setlength{\listparindent}{0cm}\setlength{\topsep}{0cm}%
      \setlength{\baselineskip}{0.5cm}\setlength{\partopsep}{0cm}%
      \setlength{\parsep}{0cm}\setlength{\parindent}{0cm}%
      \setlength{\itemsep}{1.5cm plus0.5cm minus0.5cm}%
    }
  }{\end{list}}
\newcommand{\bioplot}[2]{%
  \begin{figure*}[htbp]
    \setlength{\fboxsep}{0.5ex}
    \makebox{\parbox{0.71\textwidth}{#2}
      \hfill
      \begin{minipage}{0.29\textwidth}
        \protect\centering
        \fbox{\includegraphics[height=1.6in]{#1}}
      \end{minipage}}
  \end{figure*}}
\newcommand{\ps@georeport}{%
\renewcommand{\@oddhead}{%
\seg@foot\hfill\textit{\seg@rhead}\hfill\thepage}
\renewcommand{\@evenhead}{%
\thepage\hfill\textit{\seg@lhead}\hfill\seg@foot}
\renewcommand{\@oddfoot}{}\renewcommand{\@evenfoot}{}}
\newcommand{\ps@geobook}{%
\renewcommand{\@oddfoot}{\hfill\textit{\seg@foot}\hfill}%
\renewcommand{\@evenfoot}{\@oddfoot}}
\ifthenelse{\boolean{@book}}{%
\newcommand{\ps@segtwo}{\ps@geobook}\pagestyle{geobook}}{%
\newcommand{\ps@segtwo}{%
  \renewcommand{\@oddfoot}{\hfill\thepage\hfill}%
  \renewcommand{\@evenfoot}{\@oddfoot}
  \renewcommand{\@oddhead}{%
    \footnotesize\textit{\geo@header\seg@head}\hfill\normalsize}%
  \renewcommand{\@evenhead}{\@oddhead}}\pagestyle{georeport}}
\let\seg@large\large
\let\seg@Large\Large
\renewcommand{\Large}{\protect\seg@large}
\let\seg@section\section
\let\seg@subsection\subsection
\let\seg@subsubsection\subsubsection
\newcommand{\geosection}{%
\@startsection {section}{1}{\z@}%
{-3.5ex \@plus -1ex \@minus -.2ex}%
{1ex \@plus .2ex}%
{\centering%
\normalfont\ifthenelse{\boolean{@twoc}}{}{\Large}\bfseries}}
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
\ifthenelse{\boolean{@book}}{}{%
\def\section{\@ifstar{\segsection*}{\segsection*}}
\def\subsection{\@ifstar{\seg@subsection*}{\seg@subsection*}}%
\def\subsubsection{\@ifstar{\seg@subsubsection*}{\seg@subsubsection*}}}
\newcounter{@append}
\providecommand{\append@name}{A}
\providecommand{\appendname}[1]{\renewcommand{\append@name}{#1}}
\renewcommand{\appendix}{%
  \ifthenelse{\equal{\append@name}{A}}{%
    \renewcommand{\append@name}{\Alph{@append}}}{}
  \setcounter{equation}{0}%
  \renewcommand{\theequation}{\mbox{\append@name-\arabic{equation}}}%
  \setcounter{figure}{0}%
  \renewcommand{\thefigure}{\append@name-\arabic{figure}}%
  \stepcounter{@append}
}
\providecommand{\append}[1]{%
  \appendix\section{Appendix \append@name}\section{#1}}
\RequirePackage{natbib}
\bibpunct[,]{(}{)}{;}{a}{,}{,}
\DeclareRobustCommand{\refname}{REFERENCES}
\ifthenelse{\boolean{@book}}
{\renewcommand\bibsection{\section*{\refname}}}
{\renewcommand\bibsection{\section{\refname}}}
\RequirePackage{listings}
\renewcommand{\lst@inputpath}{./\paper@dir/}
\endinput
%%
%% End of file `georeport.cls'.
