int igetx, igety;
#define GETXY(XXX,YYY)	igetx = geth(pltin);\
			igety = geth(pltin);\
			vptodevxy(igetx,igety,&igetx,&igety);\
		        XXX = igetx;\
		        YYY = igety

#define GETXY_TEXT(XXX,YYY)	igetx = geth(pltin);\
				igety = geth(pltin);\
				vptodevxy_text(igetx,igety,&igetx,&igety);\
		        	XXX = igetx;\
		        	YYY = igety
