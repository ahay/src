\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{handout}[2003/09/01 v0.1 handouts for imaging class]

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions*

\LoadClass{article}

\RequirePackage{times}
\RequirePackage{mathptm}

\setlength{\textheight}{9.5in}
\setlength{\textwidth}{6.5in}
\setlength{\parskip}{1.5ex plus0.1ex minus0.1ex} % should be a *rubber* length
\setlength{\oddsidemargin}{.01in}
\setlength{\evensidemargin}{.01in}
\setlength{\topmargin}{-0.51in}
\setlength{\textfloatsep}{30pt plus 3pt minus 6pt}
\setlength{\headsep}{0in}

\newcommand{\@name}{Sergey Fomel}
\renewcommand{\author}[1]{\renewcommand{\@name}{#1}}

\renewcommand{\title}[1]{%
  \begin{center}
{\small
  GEO391 Seismic Wavefield Imaging \hfill #1 \\
  \@name \hfill \today} \\ 
\hrulefill
\end{center}}

% to make section headings smaller
\renewcommand{\Large}{\protect\normalsize\centering}%
\renewcommand{\large}{\protect\normalsize}

\let\sep@section\section
\let\sep@subsection\subsection
\let\sep@subsubsection\subsubsection
\newcommand{\sepsection}[1]{\sep@section*{#1}}
\renewcommand{\section}[2][]{\sepsection{\MakeUppercase{#2}}}%
\def\subsection{\sep@subsection*}%
\def\subsubsection{\sep@subsubsection*}

\RequirePackage{natbib}
\DeclareRobustCommand{\refname}{REFERENCES}
\DeclareRobustCommand{\bibname}{\protect\large~REFERENCES}
\renewcommand{\bibsection}{\section{\refname}}

\RequirePackage{graphicx}

\newcommand{\figdir}{.}
\newcommand{\inputdir}[1]{\renewcommand{\figdir}{#1}}

\newcommand{\plot}[4][htbp]{%
\begin{figure*}[#1]
\centering
\includegraphics[#3]{\figdir/Fig/#2}
\caption{#4}
\label{fig:#2}
\end{figure*}}

\newcommand{\sideplot}[4][htbp]{%
\begin{figure*}[#1]
\makebox{%
  \begin{minipage}{0.4\textwidth}
    \caption{#4}
    \label{fig:#2}
  \end{minipage}
  \hfill
  \begin{minipage}{0.6\textwidth}
    \centering
    \includegraphics[#3]{Fig/#2}
  \end{minipage}}
\end{figure*}}

\newlength{\box@width}
\setlength{\box@width}{\textwidth}
\addtolength{\box@width}{-20pt} 
\setlength{\fboxsep}{10pt}

\newcommand{\boxit}[1]{\noindent\fbox{\parbox{\box@width}{#1}}}
