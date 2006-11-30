\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{pcsslides}[2005/12/16 v0.1 slides for class]

\RequirePackage{ifthen}

%\usetheme[hideothersubsections]{PCS}
\LoadClass{beamer}

\DeclareGraphicsExtensions{.pdf,.png}

% path to figures
\newcommand{\@path}{.}
\newcommand{\inputdir}[1]{\renewcommand{\@path}{#1}}
\providecommand{\figdir}{Fig}
\providecommand{\fig@file}[1]{\@path/\figdir/#1}

% plot command
\newcommand{\plot}[4][htbp]{
  \begin{figure}[#1]
    \centering
    \includegraphics[#3]{\fig@file{#2}}
    \caption{#4}
    \label{fig:#2}
  \end{figure}
}

% sideplot command
\newcommand{\sideplot}[4][htbp]{
  \begin{figure}[#1]
    \makebox{
      \begin{minipage}{0.4\textwidth}
	\caption{#4}
	\label{fig:#2}
      \end{minipage}
      \hfill
      \begin{minipage}{0.6\textwidth}
	\centering
	\includegraphics[#3]{\fig@file{#2}}
    \end{minipage}}
\end{figure}}

% multiplot
\RequirePackage{subfigure}

\def\next@item#1,#2?{#1}
\def\rest@item#1,#2?{#2}
\newcounter{sub@fig}
\newcommand{\multiplot}[5][htbp]{
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
      \label{fig:\item@i}}
    \ifthenelse{\value{sub@fig} = #2}{\\ \setcounter{sub@fig}{0}}{}
    \repeat
    \caption{#5}
    \label{fig:#3}
  \end{figure}
}

% simple figure
\newcommand{\figplot}[4][htbp]{
  \begin{figure}[#1]
    \centering
    \includegraphics[#3]{\fig@file{#2}}
  \end{figure}
}

% plot using overpic
\newcommand{\oplot}[3]{
  \begin{center}
    \begin{overpic}[#2] {\@path/\figdir/#1}
      \begin{picture}(1,1) #3 \end{picture}
    \end{overpic}
  \end{center}
}

\newcommand{\otwoplot}[3]{
  \begin{columns}[t]
    \column{.5\textwidth} \oplot{#1}{#3}
    \column{.5\textwidth} \oplot{#2}{#3}
  \end{columns}
}

% ------------------------------------------------------------
\newcommand{\pcscover}{
  \begin{frame}
    \titlepage       % generate title page
  \end{frame}
}

\newcommand{\pcsoutline}{
  \AtBeginSection[]{
    \begin{frame}<beamer>
      \frametitle{Outline}             % add outline slide
      \tableofcontents[currentsection] % highlight current section
    \end{frame}
  }
}
% ------------------------------------------------------------
