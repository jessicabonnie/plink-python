require(tikzDevice)  # So that default options are set
options(tikzLatexPackages = c( 
  getOption('tikzLatexPackages'),  # The original contents: required stuff
  \def\tooltiptarget{\phantom{\rule{1mm}{1mm}}}

\newbox\tempboxa
\setbox\tempboxa=\hbox{} 
\immediate\pdfxform\tempboxa 
\edef\emptyicon{\the\pdflastxform}

\newcommand\tooltip[1]{%
  \pdfstartlink user{%
    /Subtype /Text
    /Contents  (#1)
    /AP <<
      /N \emptyicon\space 0 R
    >>
  }%
  \tooltiptarget%
  \pdfendlink%
}

))
