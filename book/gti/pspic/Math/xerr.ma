Approx[a_,g_]:=Sqrt[Cos[a]^2/3 + 
  (2*Sqrt[-g^2 + Cos[a]^4 + 3*g*Cos[a]*Sin[a]])/3];
Exact[a_,g_]:=-(-3*Cos[a]^2*Sin[a] - 2*Sin[a]^3)/(12*g) + 
 Sqrt[-(-48*g^2*Cos[a]^2 - 4*Cos[a]^6 - 6*g^2*Sin[a]^2 - 
       3*Cos[a]^4*Sin[a]^2)/(36*g^2) - 
    (48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 6*g^2*Sin[a]^2 + 
      3*Cos[a]^4*Sin[a]^2)/(108*g^2) + 
    (-3*Cos[a]^2*Sin[a] - 2*Sin[a]^3)^2/(36*g^2) + 
    (2^(1/3)*(3072*g^6 + 6912*g^4*Cos[a]^4 + 
       2112*g^2*Cos[a]^8 + 16*Cos[a]^12 + 
       12672*g^4*Cos[a]^2*Sin[a]^2 + 5952*g^2*Cos[a]^6*
        Sin[a]^2 + 24*Cos[a]^10*Sin[a]^2 + 
       6480*g^4*Sin[a]^4 + 5832*g^2*Cos[a]^4*Sin[a]^4 + 
       9*Cos[a]^8*Sin[a]^4 + 1944*g^2*Cos[a]^2*
        Sin[a]^6))/(3*g^2*
      (-93312*(48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
           6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)^3 - 
        5038848*(48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
          6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)*
         (3*g*Cos[a]^2*Sin[a] + 2*g*Sin[a]^3)*
         (-32*g^3*Sin[a] + 40*g*Cos[a]^4*Sin[a] + 
          27*g*Cos[a]^2*Sin[a]^3) + 45349632*g^2*
         (-32*g^3*Sin[a] + 40*g*Cos[a]^4*Sin[a] + 
           27*g*Cos[a]^2*Sin[a]^3)^2 + 3359232*g^2*
         (48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
          6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)*
         (256*g^4 + 384*g^2*Cos[a]^4 + 144*Cos[a]^8 + 
          1296*g^2*Cos[a]^2*Sin[a]^2 + 108*Cos[a]^6*
           Sin[a]^2 + 729*g^2*Sin[a]^4) + 
        5038848*(3*g*Cos[a]^2*Sin[a] + 2*g*Sin[a]^3)^2*
         (256*g^4 + 384*g^2*Cos[a]^4 + 144*Cos[a]^8 + 
          1296*g^2*Cos[a]^2*Sin[a]^2 + 108*Cos[a]^6*
           Sin[a]^2 + 729*g^2*Sin[a]^4) + 
        Sqrt[-4*(3981312*g^6 + 8957952*g^4*Cos[a]^4 + 
             2737152*g^2*Cos[a]^8 + 20736*Cos[a]^12 + 
             16422912*g^4*Cos[a]^2*Sin[a]^2 + 
             7713792*g^2*Cos[a]^6*Sin[a]^2 + 
             31104*Cos[a]^10*Sin[a]^2 + 8398080*g^4*
              Sin[a]^4 + 7558272*g^2*Cos[a]^4*Sin[a]^4 + 
             11664*Cos[a]^8*Sin[a]^4 + 2519424*g^2*
              Cos[a]^2*Sin[a]^6)^3 + 
          (-93312*(48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 6*g^2*
                Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)^3 - 
            5038848*(48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
              6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)*
             (3*g*Cos[a]^2*Sin[a] + 2*g*Sin[a]^3)*
             (-32*g^3*Sin[a] + 40*g*Cos[a]^4*Sin[a] + 
              27*g*Cos[a]^2*Sin[a]^3) + 45349632*g^2*
             (-32*g^3*Sin[a] + 40*g*Cos[a]^4*Sin[a] + 27*
                g*Cos[a]^2*Sin[a]^3)^2 + 3359232*g^2*
             (48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
              6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)*
             (256*g^4 + 384*g^2*Cos[a]^4 + 144*Cos[a]^
                8 + 1296*g^2*Cos[a]^2*Sin[a]^2 + 
              108*Cos[a]^6*Sin[a]^2 + 729*g^2*Sin[a]^
                4) + 5038848*(3*g*Cos[a]^2*Sin[a] + 2*g*
                Sin[a]^3)^2*(256*g^4 + 384*g^2*Cos[a]^
                4 + 144*Cos[a]^8 + 1296*g^2*Cos[a]^2*
               Sin[a]^2 + 108*Cos[a]^6*Sin[a]^2 + 
              729*g^2*Sin[a]^4))^2])^(1/3)) + 
    (-93312*(48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
          6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)^3 - 
       5038848*(48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
         6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)*
        (3*g*Cos[a]^2*Sin[a] + 2*g*Sin[a]^3)*
        (-32*g^3*Sin[a] + 40*g*Cos[a]^4*Sin[a] + 
         27*g*Cos[a]^2*Sin[a]^3) + 45349632*g^2*
        (-32*g^3*Sin[a] + 40*g*Cos[a]^4*Sin[a] + 
          27*g*Cos[a]^2*Sin[a]^3)^2 + 3359232*g^2*
        (48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 6*g^2*Sin[a]^2 + 
         3*Cos[a]^4*Sin[a]^2)*(256*g^4 + 
         384*g^2*Cos[a]^4 + 144*Cos[a]^8 + 
         1296*g^2*Cos[a]^2*Sin[a]^2 + 108*Cos[a]^6*
          Sin[a]^2 + 729*g^2*Sin[a]^4) + 
       5038848*(3*g*Cos[a]^2*Sin[a] + 2*g*Sin[a]^3)^2*
        (256*g^4 + 384*g^2*Cos[a]^4 + 144*Cos[a]^8 + 
         1296*g^2*Cos[a]^2*Sin[a]^2 + 108*Cos[a]^6*
          Sin[a]^2 + 729*g^2*Sin[a]^4) + 
       Sqrt[-4*(3981312*g^6 + 8957952*g^4*Cos[a]^4 + 
            2737152*g^2*Cos[a]^8 + 20736*Cos[a]^12 + 
            16422912*g^4*Cos[a]^2*Sin[a]^2 + 
            7713792*g^2*Cos[a]^6*Sin[a]^2 + 
            31104*Cos[a]^10*Sin[a]^2 + 8398080*g^4*
             Sin[a]^4 + 7558272*g^2*Cos[a]^4*Sin[a]^4 + 
            11664*Cos[a]^8*Sin[a]^4 + 2519424*g^2*
             Cos[a]^2*Sin[a]^6)^3 + 
         (-93312*(48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
              6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)^3 - 
           5038848*(48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
             6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)*
            (3*g*Cos[a]^2*Sin[a] + 2*g*Sin[a]^3)*
            (-32*g^3*Sin[a] + 40*g*Cos[a]^4*Sin[a] + 
             27*g*Cos[a]^2*Sin[a]^3) + 45349632*g^2*
            (-32*g^3*Sin[a] + 40*g*Cos[a]^4*Sin[a] + 
              27*g*Cos[a]^2*Sin[a]^3)^2 + 3359232*g^2*
            (48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
             6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)*
            (256*g^4 + 384*g^2*Cos[a]^4 + 144*Cos[a]^8 + 
             1296*g^2*Cos[a]^2*Sin[a]^2 + 108*Cos[a]^6*
              Sin[a]^2 + 729*g^2*Sin[a]^4) + 
           5038848*(3*g*Cos[a]^2*Sin[a] + 2*g*Sin[a]^3)^
             2*(256*g^4 + 384*g^2*Cos[a]^4 + 
             144*Cos[a]^8 + 1296*g^2*Cos[a]^2*Sin[a]^2 + 
             108*Cos[a]^6*Sin[a]^2 + 729*g^2*Sin[a]^4))^
          2])^(1/3)/(3888*2^(1/3)*g^2)]/2 - 
 Sqrt[-(-48*g^2*Cos[a]^2 - 4*Cos[a]^6 - 6*g^2*Sin[a]^2 - 
       3*Cos[a]^4*Sin[a]^2)/(36*g^2) + 
    (48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 6*g^2*Sin[a]^2 + 
      3*Cos[a]^4*Sin[a]^2)/(108*g^2) + 
    (-3*Cos[a]^2*Sin[a] - 2*Sin[a]^3)^2/(18*g^2) - 
    (2^(1/3)*(3072*g^6 + 6912*g^4*Cos[a]^4 + 
       2112*g^2*Cos[a]^8 + 16*Cos[a]^12 + 
       12672*g^4*Cos[a]^2*Sin[a]^2 + 5952*g^2*Cos[a]^6*
        Sin[a]^2 + 24*Cos[a]^10*Sin[a]^2 + 
       6480*g^4*Sin[a]^4 + 5832*g^2*Cos[a]^4*Sin[a]^4 + 
       9*Cos[a]^8*Sin[a]^4 + 1944*g^2*Cos[a]^2*
        Sin[a]^6))/(3*g^2*
      (-93312*(48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
           6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)^3 - 
        5038848*(48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
          6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)*
         (3*g*Cos[a]^2*Sin[a] + 2*g*Sin[a]^3)*
         (-32*g^3*Sin[a] + 40*g*Cos[a]^4*Sin[a] + 
          27*g*Cos[a]^2*Sin[a]^3) + 45349632*g^2*
         (-32*g^3*Sin[a] + 40*g*Cos[a]^4*Sin[a] + 
           27*g*Cos[a]^2*Sin[a]^3)^2 + 3359232*g^2*
         (48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
          6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)*
         (256*g^4 + 384*g^2*Cos[a]^4 + 144*Cos[a]^8 + 
          1296*g^2*Cos[a]^2*Sin[a]^2 + 108*Cos[a]^6*
           Sin[a]^2 + 729*g^2*Sin[a]^4) + 
        5038848*(3*g*Cos[a]^2*Sin[a] + 2*g*Sin[a]^3)^2*
         (256*g^4 + 384*g^2*Cos[a]^4 + 144*Cos[a]^8 + 
          1296*g^2*Cos[a]^2*Sin[a]^2 + 108*Cos[a]^6*
           Sin[a]^2 + 729*g^2*Sin[a]^4) + 
        Sqrt[-4*(3981312*g^6 + 8957952*g^4*Cos[a]^4 + 
             2737152*g^2*Cos[a]^8 + 20736*Cos[a]^12 + 
             16422912*g^4*Cos[a]^2*Sin[a]^2 + 
             7713792*g^2*Cos[a]^6*Sin[a]^2 + 
             31104*Cos[a]^10*Sin[a]^2 + 8398080*g^4*
              Sin[a]^4 + 7558272*g^2*Cos[a]^4*Sin[a]^4 + 
             11664*Cos[a]^8*Sin[a]^4 + 2519424*g^2*
              Cos[a]^2*Sin[a]^6)^3 + 
          (-93312*(48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 6*g^2*
                Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)^3 - 
            5038848*(48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
              6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)*
             (3*g*Cos[a]^2*Sin[a] + 2*g*Sin[a]^3)*
             (-32*g^3*Sin[a] + 40*g*Cos[a]^4*Sin[a] + 
              27*g*Cos[a]^2*Sin[a]^3) + 45349632*g^2*
             (-32*g^3*Sin[a] + 40*g*Cos[a]^4*Sin[a] + 27*
                g*Cos[a]^2*Sin[a]^3)^2 + 3359232*g^2*
             (48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
              6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)*
             (256*g^4 + 384*g^2*Cos[a]^4 + 144*Cos[a]^
                8 + 1296*g^2*Cos[a]^2*Sin[a]^2 + 
              108*Cos[a]^6*Sin[a]^2 + 729*g^2*Sin[a]^
                4) + 5038848*(3*g*Cos[a]^2*Sin[a] + 2*g*
                Sin[a]^3)^2*(256*g^4 + 384*g^2*Cos[a]^
                4 + 144*Cos[a]^8 + 1296*g^2*Cos[a]^2*
               Sin[a]^2 + 108*Cos[a]^6*Sin[a]^2 + 
              729*g^2*Sin[a]^4))^2])^(1/3)) - 
    (-93312*(48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
          6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)^3 - 
       5038848*(48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
         6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)*
        (3*g*Cos[a]^2*Sin[a] + 2*g*Sin[a]^3)*
        (-32*g^3*Sin[a] + 40*g*Cos[a]^4*Sin[a] + 
         27*g*Cos[a]^2*Sin[a]^3) + 45349632*g^2*
        (-32*g^3*Sin[a] + 40*g*Cos[a]^4*Sin[a] + 
          27*g*Cos[a]^2*Sin[a]^3)^2 + 3359232*g^2*
        (48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 6*g^2*Sin[a]^2 + 
         3*Cos[a]^4*Sin[a]^2)*(256*g^4 + 
         384*g^2*Cos[a]^4 + 144*Cos[a]^8 + 
         1296*g^2*Cos[a]^2*Sin[a]^2 + 108*Cos[a]^6*
          Sin[a]^2 + 729*g^2*Sin[a]^4) + 
       5038848*(3*g*Cos[a]^2*Sin[a] + 2*g*Sin[a]^3)^2*
        (256*g^4 + 384*g^2*Cos[a]^4 + 144*Cos[a]^8 + 
         1296*g^2*Cos[a]^2*Sin[a]^2 + 108*Cos[a]^6*
          Sin[a]^2 + 729*g^2*Sin[a]^4) + 
       Sqrt[-4*(3981312*g^6 + 8957952*g^4*Cos[a]^4 + 
            2737152*g^2*Cos[a]^8 + 20736*Cos[a]^12 + 
            16422912*g^4*Cos[a]^2*Sin[a]^2 + 
            7713792*g^2*Cos[a]^6*Sin[a]^2 + 
            31104*Cos[a]^10*Sin[a]^2 + 8398080*g^4*
             Sin[a]^4 + 7558272*g^2*Cos[a]^4*Sin[a]^4 + 
            11664*Cos[a]^8*Sin[a]^4 + 2519424*g^2*
             Cos[a]^2*Sin[a]^6)^3 + 
         (-93312*(48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
              6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)^3 - 
           5038848*(48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
             6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)*
            (3*g*Cos[a]^2*Sin[a] + 2*g*Sin[a]^3)*
            (-32*g^3*Sin[a] + 40*g*Cos[a]^4*Sin[a] + 
             27*g*Cos[a]^2*Sin[a]^3) + 45349632*g^2*
            (-32*g^3*Sin[a] + 40*g*Cos[a]^4*Sin[a] + 
              27*g*Cos[a]^2*Sin[a]^3)^2 + 3359232*g^2*
            (48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
             6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)*
            (256*g^4 + 384*g^2*Cos[a]^4 + 144*Cos[a]^8 + 
             1296*g^2*Cos[a]^2*Sin[a]^2 + 108*Cos[a]^6*
              Sin[a]^2 + 729*g^2*Sin[a]^4) + 
           5038848*(3*g*Cos[a]^2*Sin[a] + 2*g*Sin[a]^3)^
             2*(256*g^4 + 384*g^2*Cos[a]^4 + 
             144*Cos[a]^8 + 1296*g^2*Cos[a]^2*Sin[a]^2 + 
             108*Cos[a]^6*Sin[a]^2 + 729*g^2*Sin[a]^4))^
          2])^(1/3)/(3888*2^(1/3)*g^2) + 
    (((-48*g^2*Cos[a]^2 - 4*Cos[a]^6 - 6*g^2*Sin[a]^2 - 
         3*Cos[a]^4*Sin[a]^2)*(-3*Cos[a]^2*Sin[a] - 
         2*Sin[a]^3))/(27*g^3) - 
      (-3*Cos[a]^2*Sin[a] - 2*Sin[a]^3)^3/(27*g^3) - 
      (2*(-32*g^2*Sin[a] + 40*Cos[a]^4*Sin[a] + 
         27*Cos[a]^2*Sin[a]^3))/(9*g))/
     (4*Sqrt[-(-48*g^2*Cos[a]^2 - 4*Cos[a]^6 - 
           6*g^2*Sin[a]^2 - 3*Cos[a]^4*Sin[a]^2)/
         (36*g^2) - (48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
          6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)/
         (108*g^2) + (-3*Cos[a]^2*Sin[a] - 2*Sin[a]^3)^2/
         (36*g^2) + (2^(1/3)*(3072*g^6 + 
           6912*g^4*Cos[a]^4 + 2112*g^2*Cos[a]^8 + 
           16*Cos[a]^12 + 12672*g^4*Cos[a]^2*Sin[a]^2 + 
           5952*g^2*Cos[a]^6*Sin[a]^2 + 24*Cos[a]^10*
            Sin[a]^2 + 6480*g^4*Sin[a]^4 + 
           5832*g^2*Cos[a]^4*Sin[a]^4 + 9*Cos[a]^8*
            Sin[a]^4 + 1944*g^2*Cos[a]^2*Sin[a]^6))/
         (3*g^2*(-93312*(48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
               6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)^3 - 
            5038848*(48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
              6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)*
             (3*g*Cos[a]^2*Sin[a] + 2*g*Sin[a]^3)*
             (-32*g^3*Sin[a] + 40*g*Cos[a]^4*Sin[a] + 
              27*g*Cos[a]^2*Sin[a]^3) + 45349632*g^2*
             (-32*g^3*Sin[a] + 40*g*Cos[a]^4*Sin[a] + 27*
                g*Cos[a]^2*Sin[a]^3)^2 + 3359232*g^2*
             (48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
              6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)*
             (256*g^4 + 384*g^2*Cos[a]^4 + 144*Cos[a]^
                8 + 1296*g^2*Cos[a]^2*Sin[a]^2 + 
              108*Cos[a]^6*Sin[a]^2 + 729*g^2*Sin[a]^
                4) + 5038848*(3*g*Cos[a]^2*Sin[a] + 2*g*
                Sin[a]^3)^2*(256*g^4 + 384*g^2*Cos[a]^
                4 + 144*Cos[a]^8 + 1296*g^2*Cos[a]^2*
               Sin[a]^2 + 108*Cos[a]^6*Sin[a]^2 + 
              729*g^2*Sin[a]^4) + Sqrt[
             -4*(3981312*g^6 + 8957952*g^4*Cos[a]^4 + 
                 2737152*g^2*Cos[a]^8 + 20736*Cos[a]^
                   12 + 16422912*g^4*Cos[a]^2*Sin[a]^2 + 
                 7713792*g^2*Cos[a]^6*Sin[a]^2 + 
                 31104*Cos[a]^10*Sin[a]^2 + 8398080*g^4*
                  Sin[a]^4 + 7558272*g^2*Cos[a]^4*
                  Sin[a]^4 + 11664*Cos[a]^8*Sin[a]^4 + 
                 2519424*g^2*Cos[a]^2*Sin[a]^6)^3 + 
              (-93312*(48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
                   6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)^
                  3 - 5038848*(48*g^2*Cos[a]^2 + 
                  4*Cos[a]^6 + 6*g^2*Sin[a]^2 + 
                  3*Cos[a]^4*Sin[a]^2)*(3*g*Cos[a]^2*
                   Sin[a] + 2*g*Sin[a]^3)*(-32*g^3*
                   Sin[a] + 40*g*Cos[a]^4*Sin[a] + 
                  27*g*Cos[a]^2*Sin[a]^3) + 45349632*g^2*
                 (-32*g^3*Sin[a] + 40*g*Cos[a]^4*
                   Sin[a] + 27*g*Cos[a]^2*Sin[a]^3)^2 + 
                3359232*g^2*(48*g^2*Cos[a]^2 + 
                  4*Cos[a]^6 + 6*g^2*Sin[a]^2 + 
                  3*Cos[a]^4*Sin[a]^2)*(256*g^4 + 
                  384*g^2*Cos[a]^4 + 144*Cos[a]^8 + 
                  1296*g^2*Cos[a]^2*Sin[a]^2 + 
                  108*Cos[a]^6*Sin[a]^2 + 729*g^2*
                   Sin[a]^4) + 5038848*(3*g*Cos[a]^2*
                   Sin[a] + 2*g*Sin[a]^3)^2*(256*g^4 + 
                  384*g^2*Cos[a]^4 + 144*Cos[a]^8 + 
                  1296*g^2*Cos[a]^2*Sin[a]^2 + 
                  108*Cos[a]^6*Sin[a]^2 + 729*g^2*
                   Sin[a]^4))^2])^(1/3)) + 
        (-93312*(48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
              6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)^3 - 
           5038848*(48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
             6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)*
            (3*g*Cos[a]^2*Sin[a] + 2*g*Sin[a]^3)*
            (-32*g^3*Sin[a] + 40*g*Cos[a]^4*Sin[a] + 
             27*g*Cos[a]^2*Sin[a]^3) + 45349632*g^2*
            (-32*g^3*Sin[a] + 40*g*Cos[a]^4*Sin[a] + 
              27*g*Cos[a]^2*Sin[a]^3)^2 + 3359232*g^2*
            (48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
             6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)*
            (256*g^4 + 384*g^2*Cos[a]^4 + 144*Cos[a]^8 + 
             1296*g^2*Cos[a]^2*Sin[a]^2 + 108*Cos[a]^6*
              Sin[a]^2 + 729*g^2*Sin[a]^4) + 
           5038848*(3*g*Cos[a]^2*Sin[a] + 2*g*Sin[a]^3)^
             2*(256*g^4 + 384*g^2*Cos[a]^4 + 
             144*Cos[a]^8 + 1296*g^2*Cos[a]^2*Sin[a]^2 + 
             108*Cos[a]^6*Sin[a]^2 + 729*g^2*Sin[a]^4) + 
           Sqrt[-4*(3981312*g^6 + 8957952*g^4*Cos[a]^4 + 
                2737152*g^2*Cos[a]^8 + 20736*Cos[a]^12 + 
                16422912*g^4*Cos[a]^2*Sin[a]^2 + 
                7713792*g^2*Cos[a]^6*Sin[a]^2 + 
                31104*Cos[a]^10*Sin[a]^2 + 8398080*g^4*
                 Sin[a]^4 + 7558272*g^2*Cos[a]^4*
                 Sin[a]^4 + 11664*Cos[a]^8*Sin[a]^4 + 
                2519424*g^2*Cos[a]^2*Sin[a]^6)^3 + 
             (-93312*(48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
                  6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)^
                 3 - 5038848*(48*g^2*Cos[a]^2 + 
                 4*Cos[a]^6 + 6*g^2*Sin[a]^2 + 
                 3*Cos[a]^4*Sin[a]^2)*(3*g*Cos[a]^2*
                  Sin[a] + 2*g*Sin[a]^3)*(-32*g^3*
                  Sin[a] + 40*g*Cos[a]^4*Sin[a] + 
                 27*g*Cos[a]^2*Sin[a]^3) + 45349632*g^2*
                (-32*g^3*Sin[a] + 40*g*Cos[a]^4*Sin[a] + 
                  27*g*Cos[a]^2*Sin[a]^3)^2 + 3359232*
                g^2*(48*g^2*Cos[a]^2 + 4*Cos[a]^6 + 
                 6*g^2*Sin[a]^2 + 3*Cos[a]^4*Sin[a]^2)*
                (256*g^4 + 384*g^2*Cos[a]^4 + 
                 144*Cos[a]^8 + 1296*g^2*Cos[a]^2*
                  Sin[a]^2 + 108*Cos[a]^6*Sin[a]^2 + 
                 729*g^2*Sin[a]^4) + 5038848*
                (3*g*Cos[a]^2*Sin[a] + 2*g*Sin[a]^3)^2*
                (256*g^4 + 384*g^2*Cos[a]^4 + 
                 144*Cos[a]^8 + 1296*g^2*Cos[a]^2*
                  Sin[a]^2 + 108*Cos[a]^6*Sin[a]^2 + 
                 729*g^2*Sin[a]^4))^2])^(1/3)/
         (3888*2^(1/3)*g^2)])]/2;
err2[a_,g_]:=(Approx[a,g]-Exact[a,g])/Exact[a,g];
err1[a_,g_]:=(Cos[a]-Exact[a,g])/Exact[a,g];
errplot[p_,label_]:=Plot[{100 Abs[err2[a Pi/180, p/100]], 
		          100 Abs[err1[a Pi/180, p/100]]}, 
			  {a, 0.01,80}, 
  Frame -> True, 
  PlotStyle -> {{Thickness[0.01]}, {Thickness[0.01],Dashing[{0.03}]}}, 
  GridLines -> None, PlotRange -> All, PlotLabel->label,
  FrameLabel -> {"Angle (degrees)","Error (%)"}];
errplot[0.1,"a"];
errplot[1,"b"];
errplot[10,"c"];
Show[GraphicsArray[{{%%%},{%%},{%}}],AspectRatio->3/GoldenRatio];
Display["junk_ma.eps", %, "EPS",ImageSize->300];
