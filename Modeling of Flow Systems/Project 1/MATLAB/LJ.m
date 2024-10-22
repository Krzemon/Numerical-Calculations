function f=LJ(dxy,d,A, B);
  f=dxy*48*(A/(d*d*d*d*d*d*d*d*d*d*d*d*d*d)-...
            B/(d*d*d*d*d*d*d*d));
        if f> 500 
            f=500;
        end;
        if f< -500 
            f=-500;
        end;