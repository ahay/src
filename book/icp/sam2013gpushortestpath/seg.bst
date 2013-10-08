% BibTeX `seglike' bibliography style (29-Apr-08 version)
% Adapted from the `apalike' style, 24-jan-88 version; for BibTeX version 0.99a.
% Copyright (C) 2007, 2008 Soc. of Expl. Geophys. & Joerg Schleicher,
% all rights reserved.
% Copying of this file is allowed, provided that if you make any changes at all
% you name it something other than `seglike.bst' or 'seg'bst'.
% This restriction helps ensure that all copies are identical.
% Differences between `apalike' and `alpha' are generally heralded by a `%'.
% Differences between this style and `apalike' are generally heralded by a `%%'.
% The file btxbst.doc has the documentation for alpha.bst.
% This style works with the natbib package.
%
% History:
%   15-sep-86   (SK,OP) Original 'apalike' version,
%                       by Susan King and Oren Patashnik.
%   10-nov-86   (OP)    Truncated the sort.key$ string to the correct length
%                       in bib.sort.order to eliminate error message.
%   24-jan-88   (OP)    Updated for BibTeX version 0.99a, from alpha.bst 0.99a;
%                       apalike now sorts by author, then year, then title;
%   12-aug-04   (JS)    First adapted 'seglike' version, by Joerg Schleicher.
%                       Parts from 'seg.bst' (v. 16-apr-94) by Martin Karrenbach
%                       & from 'geophys.bst' (v. 30-mar-04) by Joerg Schleicher.
%                       'seglike' follows the SEG style according to the
%                       'Instructions to Authors', Geophysics, Jan/Feb 2004.
%   29-dec-04   (JS)    Abbr. tp for Tectonophysics added, and format.author
%                       modified to correct name format
%   09-aug-05   (JS)    Corrections of inproceedings: use of booktitle
%                       (proceedingstitle has priority), priority of meeting
%                       over type, no case change of booktitle and meeting
%                       Rules updated to follow the SEG style according to the
%                       'Instructions to Authors', Geophysics, Jan/Feb 2005.
%   03-oct-05   (JS)    Additional field segeditor in book and inbook
%   07-jan-06   (JS)    Reorganized use of meetingnr, meeting, organization,
%                       proceedingstitle as an higher-priority alternative to
%                       booktitle in inproceedings
%   13-oct-06   (JS)    Added comma after three-em dash representing identical
%                       author names
%   08-mar-07   (JS)    Period changed to comma before publisher in incollection
%                       comma added after "ed." where editor instead of author
%                       is used, resulting double comma eliminated after editor
%                       in incollection
%   04-sep-07   (JS)    Comma added after first author also for two authors;
%                       edition in book right after title; note in parentheses;
%                       series, number format in book and inbook corrected;
%                       segeditor in parentheses
%   29-apr-08   (JS)    Changes in inproceedings: address suppressed,
%                       priority of proceedingstitle over booktitle, if no
%                       proceedingstitle is present, use always booktitle,
%                       organization after proceedingstitle or booktitle;
%                       in book: segeditor replaced by bookeditor; in inbook:
%                       chapter used as paper title, title as book title; in
%                       inbook and incollection: pages after publisher, field
%                       bookeditor introduced; in book, inbook, incollection,
%                       inproceedings and proceedings: publisher before volume
%   17-jun-08   (JS)    Comma added before and in two-author references
%   17-jul-09   (JS)    Last changes: added earticle and webpage fields in
%                       accordance with 'Instructions to Authors', Geophysics, 
%                       Jan/Feb 2009. Field earticle has entries doi, url,
%                       and accessed (or lastchecked for compatibility with
%                       other solutions) in addition to all conventional
%                       article entries. Field webpage has entries author,
%                       year, title, url, accessed/lastchecked, and note.
%                       Default for accessed/lastchecked is \today.
%                       THIS VERSION DOES NOT WORK WITH BIBTEX 0.98i.

ENTRY
%%                      new entry for earticle and webpage
  { accessed
    address
    author
    booktitle
    chapter
%%                      new entry for earticle
    doi
    edition
    editor
    howpublished
    institution
    journal
    key
%%                      new entry for earticle and webpage
    lastchecked
%    month		not used in apalike  %% nor seglike.bst
    meeting		%% introduced in accordance with SEG's seg.bst	
    meetingnr		%% introduced in adaptation of geophys.bst	
    note
    number
%%			since number is suppressed in normal journal
%%			citations, we need a new entry for Geophysics with
%%			its new page number system and for magazines the
%%			issues of which start with page number 1 each
    issue
    organization
    pages
    proceedingstitle	%% only in seglike.bst
    publisher
    school
%%                      since standard bib styles don't allow for both author
%%                      and editor fields in book, inbook, and incollection
%%                      we need an additional field, bookeditor
    bookeditor
    series
    title
    type
%%                      new entry for earticle and webpage
    url
    volume
    year
  }
  {}
  { label extra.label sort.label
    dash.label} 			%% addition for seglike.bst

INTEGERS { output.state before.all mid.sentence after.sentence after.block
           set.colon} 			%% addition for seglike.bst

FUNCTION {init.state.consts}
{ #0 'before.all :=
  #1 'mid.sentence :=
  #2 'after.sentence :=
  #3 'after.block :=
  #4 'set.colon :=
}

STRINGS { s t }

FUNCTION {output.nonnull}
{ 's :=
  output.state mid.sentence =
    { ", " * write$ }
    { output.state after.block =
	{ add.period$ write$
	  newline$
	  "\newblock " write$
	}
	{ output.state set.colon =
          { ": " * write$ }
	  { output.state before.all =
	      'write$
	      { add.period$ " " * write$ }
	      if$
	    }
	  if$
	}
      if$
      mid.sentence 'output.state :=
    }
  if$
  s
}

FUNCTION {output}
{ duplicate$ empty$
    'pop$
    'output.nonnull
  if$
}

FUNCTION {output.check}
{ 't :=
  duplicate$ empty$
    { pop$ "empty " t * " in " * cite$ * warning$ }
    'output.nonnull
  if$
}

%					apalike needs this function because
%					the year has special punctuation;
%					apalike ignores the month
FUNCTION {output.year.check}
{ year empty$
    { "empty year in " cite$ * warning$ }
    { write$
%%      " (" year * extra.label * ")" *  
%%                                      changed to SEG style (no parentheses)
       " " year * extra.label *

      mid.sentence 'output.state :=
    }
  if$
}

FUNCTION {output.bibitem}
{ newline$
  "\bibitem[" write$
  label write$
  "]{" write$
  cite$ write$
  "}" write$
  newline$
  ""
  before.all 'output.state :=
}

FUNCTION {fin.entry}
{ add.period$
  write$
  newline$
}

FUNCTION {new.colon}
{ output.state before.all =
    'skip$
    { set.colon 'output.state := }
  if$
}

FUNCTION {new.block}
{ output.state before.all =
    'skip$
    { after.block 'output.state := }
  if$
}

FUNCTION {new.sentence}
{ output.state after.block =
    'skip$
    { output.state before.all =
	'skip$
	{ after.sentence 'output.state := }
      if$
    }
  if$
}

FUNCTION {not}
{   { #0 }
    { #1 }
  if$
}

FUNCTION {and}
{   'skip$
    { pop$ #0 }
  if$
}

FUNCTION {or}
{   { pop$ #1 }
    'skip$
  if$
}

FUNCTION {new.block.checkb}
{ empty$
  swap$ empty$
  and
    'skip$
    'new.block
  if$
}

FUNCTION {field.or.null}
{ duplicate$ empty$
    { pop$ "" }
    'skip$
  if$
}

FUNCTION {emphasize}
{ duplicate$ empty$
    { pop$ "" }
%%    { "{\em " swap$ * "}" * }			%% emphasize is bold
    { "{\bf " swap$ * "}" * }			%% in seglike.bst
  if$
}

INTEGERS { nameptr namesleft numnames }

FUNCTION {format.names.old}
{ 's :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { s nameptr "{vv~}{ll}{, jj}{, f.}" format.name$ 't :=   % last name first
      nameptr #1 >
	{ namesleft #1 >
	    { ", " * t * }
 	   {
%%		numnames #2 >             %% always add comma after first author
		{ "," * }
%%		'skip$
%%	      if$
	      t "others" =
		{ " et~al." * }
		{ " and " * t * }
	      if$
	    }
	  if$
	}
	't
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
}

FUNCTION {format.names}
{ 's :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { numnames namesleft =
        { s nameptr "{vv~}{ll}{, jj}{, f.}" format.name$ 't := }  % last name first 
        { s nameptr "{f. }{vv~}{ll}{, jj}" format.name$ 't := }  % first name first
      if$
      nameptr #1 >
	{ namesleft #1 >
	    { ", " * t * }
	    { numnames #1 >
		{ "," * }
		'skip$
	      if$
	      t "others" =
		{ " et~al." * }
		{ " and " * t * }
	      if$
	    }
	  if$
	}
	't
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
}

FUNCTION {format.authors}
{ author empty$
    { "" }
    { dash.label "-"  =
        { "--------," }
        { author format.names
         ", " *
        }
      if$
    }
  if$
}

FUNCTION {format.key}			% this function is just for apalike
{ empty$
    { key field.or.null }
    { "" }
  if$
}

FUNCTION {format.bookeditors.book}       %% this function is just for seglike
{ bookeditor empty$                      %% it formats the bookeditor field
    { "" }                              %% allowed in book and inbook
    { " (" bookeditor format.names *
      bookeditor num.names$ #1 >
%%	{ ", editors" * }   
%%	{ ", editor" * }
        { ", eds." * }                  %% SEG style
        { ", ed." * }                   %% SEG style
      if$
      ")" *
      before.all 'output.state :=
    }
  if$
}

FUNCTION {format.editors}
{ editor empty$
    { "" }
    { editor format.names
      editor num.names$ #1 >
%%	{ ", editors" * }   
%%	{ ", editor" * }
        { ", eds.," * }                  %% SEG style
        { ", ed.," * }                   %% SEG style
      if$
    }
  if$
}

FUNCTION {format.bookeditors}
{ bookeditor empty$
    { "" }
    { bookeditor format.names
      bookeditor num.names$ #1 >
%%	{ ", editors" * }   
%%	{ ", editor" * }
        { ", eds.," * }                  %% SEG style
        { ", ed.," * }                   %% SEG style
      if$
    }
  if$
}

FUNCTION {format.chapter}
{ chapter empty$
    { "" }
    { chapter "t" change.case$ }
  if$
}

FUNCTION {format.url.access}
{ url empty$
    { "url empty in " cite$ * warning$ }
    { "\verb+" url * "+" * 
      lastchecked empty$
        { accessed empty$
            { "" "both lastchecked and accessed empty in " cite$ * warning$ }
            { ", accessed " * accessed * }
          if$
        }
        { ", accessed " * lastchecked * }
      if$
    }
  if$
}

FUNCTION {format.doi.access}
{ doi empty$
    {  url empty$
        { "" "both doi and url empty in " cite$ * warning$ }
        { "empty doi in " cite$ * warning$ }
      if$
    }
    { }
  if$
  lastchecked empty$
    { accessed empty$
        { "\today" "both lastchecked and accessed empty in " cite$ * warning$ }
        { "accessed " accessed * }
      if$
    }
    { "accessed " lastchecked * }
  if$
  url empty$
    { }
    { "; \verb+" * url * "+" * }
  if$
  doi empty$
    { }
    {  "; \verb+http://dx.doi.org/" * doi * "+" * }
  if$
}

FUNCTION {format.title}
{ title empty$
    { "" }
    { title "t" change.case$ }
  if$
}

FUNCTION {n.dashify}
{ 't :=
  ""
    { t empty$ not }
    { t #1 #1 substring$ "-" =
	{ t #1 #2 substring$ "--" = not
	    { "--" *
	      t #2 global.max$ substring$ 't :=
	    }
	    {   { t #1 #1 substring$ "-" = }
		{ "-" *
		  t #2 global.max$ substring$ 't :=
		}
	      while$
	    }
	  if$
	}
	{ t #1 #1 substring$ *
	  t #2 global.max$ substring$ 't :=
	}
      if$
    }
  while$
}

FUNCTION {format.btitle}
{ title emphasize
}

FUNCTION {tie.or.space.connect}
{ duplicate$ text.length$ #3 <
    { "~" }
    { " " }
  if$
  swap$ * *
}

FUNCTION {either.or.check}
{ empty$
    'pop$
    { "can't use both " swap$ * " fields in " * cite$ * warning$ }
  if$
}

FUNCTION {format.bvolume}
{ volume empty$
    { "" }
    { series empty$
        { ", " * volume emphasize * }
        { "volume" volume emphasize tie.or.space.connect
	  " {\it of} " * series *
        }
      if$
      "volume and number" number either.or.check
    }
  if$
}

FUNCTION {format.number.series}
{ volume empty$
    { number empty$
	{ series field.or.null }
 	{
%	  output.state mid.sentence =
%	    { "number" }
%	    { "Number" }
%	  if$
%	  number tie.or.space.connect
	  series empty$
	    { "there's a number but no series in " cite$ * warning$ }
%	    { " {\it in} " * series * }
%       Replaced by
	    { series ", No.~" * number * }
	  if$
	}
      if$
    }
    { "" }
  if$
}

FUNCTION {format.edition}
{ edition empty$
    { "" }
    { output.state mid.sentence =
	{ edition "l" change.case$ " ed." * }
	{ edition "t" change.case$ " ed." * }
      if$
    }
  if$
}

INTEGERS { multiresult }

FUNCTION {multi.page.check}
{ 't :=
  #0 'multiresult :=
    { multiresult not
      t empty$ not
      and
    }
    { t #1 #1 substring$
      duplicate$ "-" =
      swap$ duplicate$ "," =
      swap$ "+" =
      or or
	{ #1 'multiresult := }
	{ t #2 global.max$ substring$ 't := }
      if$
    }
  while$
  multiresult
}

FUNCTION {format.pages}
{ pages empty$
    { "" }
    { pages multi.page.check
	{ "" pages n.dashify tie.or.space.connect }
	{ "" pages tie.or.space.connect }
      if$
    }
  if$
}

%%					format.vol.num.pages from seg.bst
%%					number changed to issue
FUNCTION {format.vol.num.pages}
{ volume field.or.null emphasize
  issue empty$
    'skip$
    { ", no. " issue *  *
      volume empty$
        { "there's a number but no volume in " cite$ * warning$ }
        'skip$
      if$
    }
  if$
  pages empty$
    'skip$
    { duplicate$ empty$
        { pop$ format.pages }
        { ", " * pages n.dashify * }
      if$
    }
  if$
}

FUNCTION {format.chapter.pages}
{ chapter empty$
    'format.pages
    { type empty$
	{ "chapter" }
	{ type "l" change.case$ }
      if$
      chapter tie.or.space.connect
      pages empty$
	'skip$
	{ ", " * format.pages * }
      if$
    }
  if$
}

FUNCTION {format.in.ed.booktitle}
{ booktitle empty$
    { "" }
    { bookeditor empty$
%	{ "In " booktitle emphasize * }
%	{ "In " format.editors * ", " * booktitle emphasize * }
%%			%% modified according to seg.bst
        { "{\it in} " booktitle * }
        { "{\it in} " format.bookeditors * " " * booktitle * }
      if$
    }
  if$
}

FUNCTION {format.in.ed.title}                   %% used in inbook
{ title empty$
    { "" }
    { bookeditor empty$
%	{ "In " booktitle emphasize * }
%	{ "In " format.editors * ", " * booktitle emphasize * }
%%			%% modified according to seg.bst
        { "{\it in} " title * }
        { "{\it in} " format.bookeditors * " " * title * }
      if$
    }
  if$
}

FUNCTION {format.thesis.type}
{ type empty$
    'skip$
    { pop$
      type "t" change.case$
    }
  if$
}

FUNCTION {format.tr.number}
{ type empty$
    { "Technical Report" }
    'type
  if$
  number empty$
    { "t" change.case$ }
    { number tie.or.space.connect }
  if$
}

FUNCTION {format.article.crossref}
{ "In"							% this is for apalike
  " \cite{" * crossref * "}" *
}

FUNCTION {format.book.crossref}
{ volume empty$
    { "empty volume in " cite$ * "'s crossref of " * crossref * warning$
      "In "
    }
    { "Volume" volume tie.or.space.connect
      " of " *
    }
  if$
  "\cite{" * crossref * "}" *				% this is for apalike
}

FUNCTION {format.incoll.inproc.crossref}
{ "In"							% this is for apalike
  " \cite{" * crossref * "}" *
}

%%						%% from geophys.bst
INTEGERS { index }

FUNCTION {number.to.counter}
{ 't :=
  t text.length$ 'index :=
  t index #1 substring$ "1" =
    { t "st " * }
    { t index #1 substring$ "2" =
      { t "nd " * }
      { t index #1 substring$ "3" =
        { t "rd " * }
        { t "th " * }
      if$
      }
    if$
    }
  if$
}

%%						%% from geophys.bst
FUNCTION {format.booktitle.pr}                  %% modified to give booktitle
{ booktitle empty$                              %% priority over type
    'skip$
    { pop$
      booktitle                                 %% and no change of case
    }
  if$
}

%%						%% from geophys.bst
FUNCTION {format.booktitle}
{ booktitle empty$
    { "" }
    { booktitle "t" change.case$ }
  if$
}

%%						%% Put note in parentheses
FUNCTION {format.note}
{ note empty$
    { "" }
    { "(" note * ")" * }
  if$
}
%%                                                  new in seglike.bst
FUNCTION {format.meeting.seg}
{ meeting empty$                                 %% if no meeting given
    { "Annual International Meeting" }           %% "Ann. Int. Meet."
    { meeting                                    %% use meeting
    }
  if$
  organization empty$
    { "" * }
    { ", " *
      organization *
    }
  if$
  proceedingstitle empty$                        %% if no proceedingstitle
    { booktitle empty$                           %%   if no booktitle
        { pages empty$                           %%     if no pages
            { "" * }                             %%       nothing
            { ", Expanded Abstracts" * }         %%     else "Exp. Abst"
          if$
        }
        { ", " *                                 %%   else
          booktitle *                            %%     use booktitle
        }
      if$
    }
    { ", " *                                     %% else
      proceedingstitle *                         %%   use proceedingstitle
    }
  if$
}

FUNCTION {format.meeting.ad.ab}
{ pages empty$
    { "Presented at the " }
    { "" }
  if$
  meetingnr empty$                                     %% if no meetingnr given
    { meeting empty$                                   %% and no meeting given
        { proceedingstitle empty$
            { booktitle empty$
                { "empty booktitle in " cite$ * warning$
                  "" *
                }
                { booktitle * }                           %% use booktitle
              if$
            }
            { proceedingstitle * }                        %% use proc.title
          if$
          organization empty$
            { "" * }
            { ", " *
              organization *
            }
          if$
        }
        { format.meeting.seg * }                       %% else use meeting
      if$                                              %% organization,
    }                                                  %% proceedingstitle, and
    { meetingnr number.to.counter *                    %% meetingnr (if exists)
      format.meeting.seg *
    }
  if$
  publisher empty$                                     %% no publisher given 
    { "" * }                                           %%    do nothing
    { ", " *                                           %% else
    publisher *                                        %%    use publisher
  }
  if$
  pages empty$
    { "" * }
    { "," *
      format.pages *
    }
  if$
}

FUNCTION {article}
{ output.bibitem
  format.authors "author" output.check
  author format.key output				% special for
  output.year.check					% apalike
% new.block						%% no period
  format.title "title" output.check
% new.block						%% no period
  new.colon						%% colon instead
  crossref missing$
%    { journal emphasize "journal" output.check		%% no italics
    { journal "journal" output.check			%% for seglike
      format.vol.num.pages output
    }
    { format.article.crossref output.nonnull
      format.pages output
    }
  if$
  new.block
  format.note output
  fin.entry
}

FUNCTION {earticle}
{ output.bibitem
  format.authors "author" output.check
  author format.key output				% special for
  output.year.check					% apalike
% new.block						%% no period
  format.title "title" output.check
% new.block						%% no period
  new.colon						%% colon instead
  crossref missing$
%    { journal emphasize "journal" output.check		%% no italics
    { journal "journal" output.check			%% for seglike
      format.vol.num.pages output
    }
    { format.article.crossref output.nonnull
      format.pages output
    }
  if$
  format.doi.access output
  new.block
  format.note output
  fin.entry
}

FUNCTION {book}
{ output.bibitem
  author empty$
    { format.editors  "author and editor" output.check
      editor format.key output
    }
    { format.authors output.nonnull
      crossref missing$                                  %% Don't allow both
        { "author and editor" editor either.or.check }   %% author and editor
        'skip$
      if$
    }
  if$
  output.year.check                             % special for apalike
% new.block                                     %% no period for seglike
%  format.btitle "title" output.check	        %% treat book titles like 
  crossref missing$
    { format.title "title" output.check         %%  any titles in seglike
      format.edition output                     %% edition right after title
      new.colon					%% colon for seglike
      publisher "publisher" output.check
      format.bvolume output
      new.block
      format.number.series output
      format.bookeditors.book output             %% Allow bookeditor
%      new.sentence
%      address output				%% address suppressed
    }
    { new.block
      format.book.crossref output.nonnull
    }
  if$
  new.block
  format.note output
  fin.entry
}

FUNCTION {booklet}
{ output.bibitem
  format.authors output
  author format.key output				% special for
  output.year.check					% apalike
% new.block						%% no period
  format.title "title" output.check
% new.block						%% no period
  new.colon						%% colon instead
  howpublished output
%      address output				%% address suppressed
  new.block
  format.note output
  fin.entry
}

FUNCTION {inbook}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check
      editor format.key output
    }
    { format.authors output.nonnull
      crossref missing$                                  %% Don't allow both
        { "author and editor" editor either.or.check }   %% author and editor
        'skip$
      if$
    }
  if$
  output.year.check				% special for apalike
% new.block                             	%% no period
%  format.btitle "btitle" output.check          %% no special treatment,
  format.chapter "chapter" output.check	        %% treat chapter as title
  crossref missing$
    { format.in.ed.title "title" output.check   %% treat title as book title
      format.edition output
      new.colon
      publisher "publisher" output.check
      format.bvolume output
      format.number.series output
%      address output				%% address suppressed
    }
    { new.block
      format.book.crossref output.nonnull
    }
  if$
  format.pages output
  new.block
  format.note output
  fin.entry
}

FUNCTION {incollection}
{ output.bibitem
  format.authors "author" output.check
  author format.key output				% special for
  output.year.check					% apalike
% new.block						%% no period
  format.title "title" output.check
% new.block						%% no period
  crossref missing$
    { format.in.ed.booktitle "booktitle" output.check
      format.edition output
      new.colon
      publisher "publisher" output.check
      format.bvolume output
      format.number.series output
      format.chapter output
      format.pages output
%      address output				%% address suppressed
    }
    { format.incoll.inproc.crossref output.nonnull
      format.chapter.pages output
    }
  if$
  new.block
  format.note output
  fin.entry
}

FUNCTION {inproceedings}
{ output.bibitem
  format.authors "author" output.check
  author format.key output				% special for
  output.year.check					% apalike
% new.block						%% no period
  format.title "title" output.check
% new.block						%% no period
  new.colon						%% colon instead
  crossref missing$
    { format.meeting.ad.ab output			%% from geophys.bst
%%				 volume, number, series, and address
%%				 suppressed as well as order changed
%%				 all treated in format.meeting.ad.ab
%%      format.in.ed.booktitle "booktitle" output.check
%%      format.bvolume output
%%      format.number.series output
%%      format.pages output
%%      address output					% for apalike
%%      new.sentence					% there's no year
%%      organization output				% here so things
%%      publisher output				% are simpler
    }
    { format.incoll.inproc.crossref output.nonnull
      format.pages output
    }
  if$
  new.block
  format.note output
  fin.entry
}

FUNCTION {conference} { inproceedings }

FUNCTION {manual}
{ output.bibitem
  format.authors output
  author format.key output				% special for
  output.year.check					% apalike
% new.block						%% no period
%  format.btitle "title" output.check		%% treat book titles like 
  format.title "title" output.check		%%  any titles in seglike
  new.colon						%% colon instead
  organization address new.block.checkb
  organization output
  address output
  format.edition output
  new.block
  format.note output
  fin.entry
}

FUNCTION {mastersthesis}
{ output.bibitem
  format.authors "author" output.check
  author format.key output				% special for
  output.year.check					% apalike
% new.block						%% no period
  format.title "title" output.check
% new.block						%% no period
  new.colon						%% colon instead
  "Master's thesis" format.thesis.type output.nonnull
  school "school" output.check
  address output
  new.block
  format.note output
  fin.entry
}

FUNCTION {misc}
{ output.bibitem
  format.authors output
  author format.key output				% special for
  output.year.check					% apalike
% new.block						%% no period
  format.title output
% new.block						%% no period
  new.colon						%% colon instead
  howpublished output
  new.block
  format.note output
  fin.entry
}

FUNCTION {webpage}                                      %% new in seglike
{ output.bibitem
  format.authors output
  author format.key output
  output.year.check
  format.title output
  format.url.access output
  new.block
  format.note output
  fin.entry
}

FUNCTION {phdthesis}
{ output.bibitem
  format.authors "author" output.check
  author format.key output				% special for
  output.year.check					% apalike
% new.block						%% no period
%  format.btitle "title" output.check		%% treat book titles like 
  format.title "title" output.check		%% any titles in seglike
% new.block						%% no period
  new.colon						%% colon instead
  "PhD thesis" format.thesis.type output.nonnull
  school "school" output.check
  address output
  new.block
  format.note output
  fin.entry
}

FUNCTION {proceedings}
{ output.bibitem
  format.editors output
  editor format.key output				% special for
  output.year.check					% apalike
% new.block						%% no period
%  format.btitle "title" output.check		%% treat book titles like 
  format.title "title" output.check		%% any titles in seglike
%                                               %% no addres in seglike
% address output				% for apalike
  organization output				% a nonempty organization
  new.colon					% we always output
  publisher output				% here
  format.bvolume output
  format.number.series output
  new.block
  format.note output
  fin.entry
}

FUNCTION {techreport}
{ output.bibitem
  format.authors "author" output.check
  author format.key output				% special for
  output.year.check					% apalike
% new.block						%% no period
  format.title "title" output.check
% new.block						%% no period
  new.colon						%% colon instead
  format.tr.number output.nonnull
  institution "institution" output.check
  address output
  new.block
  format.note output
  fin.entry
}

FUNCTION {unpublished}
{ output.bibitem
  format.authors "author" output.check
  author format.key output				% special for
  output.year.check					% apalike
% new.block						%% no period
  format.title "title" output.check
  new.block
  note "note" output.check
  fin.entry
}

FUNCTION {default.type} { misc }

MACRO {jan} {"January"}

MACRO {feb} {"February"}

MACRO {mar} {"March"}

MACRO {apr} {"April"}

MACRO {may} {"May"}

MACRO {jun} {"June"}

MACRO {jul} {"July"}

MACRO {aug} {"August"}

MACRO {sep} {"September"}

MACRO {oct} {"October"}

MACRO {nov} {"November"}

MACRO {dec} {"December"}

%%				Journal names from geophys.bst added
MACRO {acmcs} {"ACM Computing Surveys"}

MACRO {acta} {"Acta Informatica"}

MACRO {bssa} {"Bulletin of the Seismological Society of America"}

MACRO {bstj} {"Bell Systems Technical Journal"}

MACRO {cacm} {"Communications of the ACM"}

MACRO {fb} {"First Break"}

MACRO {geo} {"Geophysics"}

MACRO {gj} {"Geophysical Journal"}

MACRO {gji} {"Geophysical Journal International"}

MACRO {gjras} {"Geophysical Journal of the Royal Astronomical Society"}

MACRO {gp} {"Geophysical Prospecting"}

MACRO {ibmjrd} {"IBM Journal of Research and Development"}

MACRO {ibmsj} {"IBM Systems Journal"}

MACRO {ieeese} {"IEEE Transactions on Software Engineering"}

MACRO {ieeetap} {"IEEE Transactions on Antennas and Propagation"}

MACRO {ieeetc} {"IEEE Transactions on Computers"}

MACRO {ieeetcad}
 {"IEEE Transactions on Computer-Aided Design of Integrated Circuits"}

MACRO {ip} {"Inverse Problems"}

MACRO {ipl} {"Information Processing Letters"}

MACRO {jacm} {"Journal of the ACM"}

MACRO {jag} {"Journal of Applied Geophysics"}

MACRO {jam} {"SIAM Journal of Applied Mathematics"}

MACRO {jap} {"Journal of Applied Physics"}

MACRO {jasa} {"Journal of the Acoustical Society of America"}

MACRO {jasaa} {"Journal of the Acoustical Society of America A"}

MACRO {jasab} {"Journal of the Acoustical Society of America B"}

MACRO {jcss} {"Journal of Computer and System Sciences"}

MACRO {jg} {"Journal of Geophysics"}

MACRO {jge} {"Journal of Geophysics and Engineering"}

MACRO {jgr} {"Journal of Geophysical Research"}

MACRO {jmp} {"Journal of Mathematical Physics"}

MACRO {josa} {"Journal of the Optical Society of America"}

MACRO {josaa} {"Journal of the Optical Society of America A"}

MACRO {josab} {"Journal of the Optical Society of America B"}

MACRO {jse} {"Journal of Seismic Exploration"}

MACRO {ogj} {"Oil and Gas Journal"}

MACRO {pag} {"Pure and Applied Geophysics"}

MACRO {pepi} {"Physics of the Earth and Planetory Interior"}

MACRO {pieee} {"Proceedings of the IEEE"}

MACRO {rgg} {"Russian Geology and Geophysics"}

MACRO {rs} {"Radio Science"}

MACRO {sd} {"Scientific Drilling"}

MACRO {scp} {"Science of Computer Programming"}

MACRO {sicomp} {"SIAM Journal on Computing"}

MACRO {sg} {"Surveys in Geophysics"}

MACRO {sgg} {"Soviet Geology and Geophysics"}

MACRO {sgeg} {"Studia Geophysica et Geodaetica"}

MACRO {tle} {"The Leading Edge"}

MACRO {tp} {"Tectonophysics"}

MACRO {tocs} {"ACM Transactions on Computer Systems"}

MACRO {tods} {"ACM Transactions on Database Systems"}

MACRO {tog} {"ACM Transactions on Graphics"}

MACRO {toms} {"ACM Transactions on Mathematical Software"}

MACRO {toois} {"ACM Transactions on Office Information Systems"}

MACRO {toplas} {"ACM Transactions on Programming Languages and Systems"}

MACRO {tcs} {"Theoretical Computer Science"}

MACRO {wm} {"Wave Motion"}

MACRO {dgg} {"Deutsche Geophysikalische Gesellschaft"}

MACRO {eage} {"EAGE"}

MACRO {eaeg} {"EAEG"}

MACRO {seg} {"SEG"}

MACRO {sbgf} {"Brazilian Geophysical Society"}

READ

FUNCTION {sortify}
{ purify$
  "l" change.case$
}

INTEGERS { len }

FUNCTION {chop.word}
{ 's :=
  'len :=
  s #1 len substring$ =
    { s len #1 + global.max$ substring$ }
    's
  if$
}

%			There are three apalike cases: one person (Jones),
%			two (Jones and de~Bruijn), and more (Jones et~al.).
%			This function is much like format.crossref.editors.
%
FUNCTION {format.lab.names}
{ 's :=
  s #1 "{vv~}{ll}" format.name$
  s num.names$ duplicate$
  #2 >
    { pop$ " et~al." * }
    { #2 <
	'skip$
	{ s #2 "{ff }{vv }{ll}{ jj}" format.name$ "others" =
	    { " et~al." * }
	    { " and " * s #2 "{vv~}{ll}" format.name$ * }
	  if$
	}
      if$
    }
  if$
}

FUNCTION {author.key.label}
{ author empty$
    { key empty$
	{ cite$ #1 #3 substring$ }
	'key					% apalike uses the whole key
      if$
    }
    { author format.lab.names }
  if$
}

FUNCTION {author.editor.key.label}
{ author empty$
    { editor empty$
	{ key empty$
	    { cite$ #1 #3 substring$ }
	    'key				% apalike uses the whole key
	  if$
	}
	{ editor format.lab.names }
      if$
    }
    { author format.lab.names }
  if$
}

FUNCTION {editor.key.label}
{ editor empty$
    { key empty$
	{ cite$ #1 #3 substring$ }
	'key			% apalike uses the whole key, no organization
      if$
    }
    { editor format.lab.names }
  if$
}

FUNCTION {calc.label}
{ type$ "book" =
  type$ "inbook" =
  or
    'author.editor.key.label
    { type$ "proceedings" =
	'editor.key.label			% apalike ignores organization
	'author.key.label			% for labeling and sorting
      if$
    }
  if$
  ", "							% these three lines are
  *							% for apalike, which
  year field.or.null purify$ #-1 #4 substring$		% uses all four digits
  *
  'label :=
}

FUNCTION {sort.format.names}
{ 's :=
  #1 'nameptr :=
  ""
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { nameptr #1 >
	{ "   " * }
	'skip$
      if$						% apalike uses initials
      s nameptr "{vv{ } }{ll{ }}{  f{ }}{  jj{ }}" format.name$ 't := % <= here
      nameptr numnames = t "others" = and
	{ "et al" * }
	{ t sortify * }
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
}

FUNCTION {sort.format.title}
{ 't :=
  "A " #2
    "An " #3
      "The " #4 t chop.word
    chop.word
  chop.word
  sortify
  #1 global.max$ substring$
}

FUNCTION {author.sort}
{ author empty$
    { key empty$
	{ "to sort, need author or key in " cite$ * warning$
	  ""
	}
	{ key sortify }
      if$
    }
    { author sort.format.names }
  if$
}

FUNCTION {author.editor.sort}
{ author empty$
    { editor empty$
	{ key empty$
	    { "to sort, need author, editor, or key in " cite$ * warning$
	      ""
	    }
	    { key sortify }
	  if$
	}
	{ editor sort.format.names }
      if$
    }
    { author sort.format.names }
  if$
}

FUNCTION {editor.sort}
{ editor empty$
    { key empty$
	{ "to sort, need editor or key in " cite$ * warning$
	  ""
	}
	{ key sortify }
      if$
    }
    { editor sort.format.names }
  if$
}

%			apalike uses two sorting passes; the first one sets the
%			labels so that the `a's, `b's, etc. can be computed;
%			the second pass puts the references in "correct" order.
%			The presort function is for the first pass. It computes
%			label, sort.label, and title, and then concatenates.
FUNCTION {presort}
{ calc.label
  label sortify
  "    "
  *
  type$ "book" =
  type$ "inbook" =
  or
    'author.editor.sort
    { type$ "proceedings" =
	'editor.sort
	'author.sort
      if$
    }
  if$
  #1 entry.max$ substring$	% for
  'sort.label :=		% apalike
  sort.label			% style
  *
  "    "
  *
  title field.or.null
  sort.format.title
  *
  #1 entry.max$ substring$
  'sort.key$ :=
}

ITERATE {presort}

SORT		% by label, sort.label, title---for final label calculation

%STRINGS { last.label next.extra }	% apalike labels are only for the text;
STRINGS { last.label next.extra
          last.sort.label}		%% addition for seglike.bst

INTEGERS { last.extra.num }		% there are none in the bibliography

FUNCTION {initialize.extra.label.stuff}	% and hence there is no `longest.label'
{ #0 int.to.chr$ 'last.label :=
  "" 'next.extra :=
  #0 'last.extra.num :=
}

FUNCTION {forward.pass}
{ last.label label =
    { last.extra.num #1 + 'last.extra.num :=
      last.extra.num int.to.chr$ 'extra.label :=
    }
    { "a" chr.to.int$ 'last.extra.num :=
      "" 'extra.label :=
      label 'last.label :=
    }
  if$
}

FUNCTION {reverse.pass}
{ next.extra "b" =
    { "a" 'extra.label := }
    'skip$
  if$
  label extra.label * 'label :=
  extra.label 'next.extra :=
}

EXECUTE {initialize.extra.label.stuff}

ITERATE {forward.pass}

REVERSE {reverse.pass}

%				Now that the label is right we sort for real,
%				on sort.label then year then title.  This is
%				for the second sorting pass.
FUNCTION {bib.sort.order}
{ sort.label
  "    "
  *
  year field.or.null sortify
  *
  "    "
  *
  title field.or.null
  sort.format.title
  *
  #1 entry.max$ substring$
  'sort.key$ :=
}

ITERATE {bib.sort.order}

SORT		% by sort.label, year, title---giving final bibliography order

%%				Now that the final order is established, we
%%				substitute multiple occurrences of the
%%				same authors by a horizontal line.
%%				(only for seglike.bst)
FUNCTION {final.pass}
{ last.sort.label sort.label =
    { "-" 'dash.label :=
    }
    { sort.label 'last.sort.label :=
      "" 'dash.label :=
    }
  if$
}

ITERATE {final.pass}

FUNCTION {begin.bib}
{ preamble$ empty$				% no \etalchar in apalike
    'skip$
    { preamble$ write$ newline$ }
  if$
  "\begin{thebibliography}{}" write$ newline$	% no labels in apalike
%%						%% normal line space
  "\itemsep0pt" write$ newline$ 		%% between references
}

EXECUTE {begin.bib}

EXECUTE {init.state.consts}

ITERATE {call.type$}

FUNCTION {end.bib}
{ newline$
  "\end{thebibliography}" write$ newline$
}

EXECUTE {end.bib}
