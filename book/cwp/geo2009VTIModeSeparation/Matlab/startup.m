

% First the fontsize.

  set(0,'defaultaxesfontsize',20);
  set(0,'defaulttextfontsize',20);

% The default figure size.  All MATLAB figures will open with this
% size and position, unless otherwise told.  This includes the Display
% Windows in dfield and pplane.

  rrs = get(0,'screensize');
  rrf = get(0,'defaultfigurepos');
  m = min((rrs(3)-200)/rrf(3), (rrs(4)-200)/rrf(4));
  defpos = [200, rrs(4)-m*rrf(4), m*rrf(3), m*rrf(4)];
  set(0,'defaultfigurepos',defpos);

% We want the line widths to be 2 points instead of the usual 1.

  set(0,'defaultaxeslinewidth',2);
  set(0,'defaultlinelinewidth',2);

% Finally define and set the global variables that provide the default
% postions for the dfield and pplane windows.

  global DFSETPOS DFKBDPOS DFSETTINGSPOS DFZOOMBACKPOS

  DFSETPOS = [15,rrs(4)-380,845,380];
  DFKBDPOS = [15,DFSETPOS(2)-240,490,230];
  DFSETTINGSPOS = [40,min(rrs(4)-420,250),820,420];
  DFZOOMBACKPOS = [329,166,717,515];


  global PPSETPOS  PPKBDPOS PPEQPTPOS PPZOOMBACKPOS
  
  PPSETPOS = [15,267,944,566];
  PPKBDPOS = [39,487,485,297];
  PPEQPTPOS = [51,315,551,395];
  PPZOOMBACKPOS = [329,166,717,515];

  global PPSETTINGSPOS PPXYPOS

  PPSETTINGSPOS = [61,135,793,375];
  PPXYPOS = [128,18,936,643];

