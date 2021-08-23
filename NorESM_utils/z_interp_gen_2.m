function [int_field] = z_interp_gen(target_bnd,field,z_bnd)

int_field=zeros(1,size(target_bnd,1));
kix=1;
for k=1:size(target_bnd,1)
 topz=target_bnd(k,1);
 botz=target_bnd(k,2);
 z_thick=0.;
 for k1=kix:size(z_bnd,1)
  if field(k1)>=0.
  if z_bnd(k1,1)>=topz & z_bnd(k1,2)<=botz & z_bnd(k1,1)<botz & z_bnd(k1,2)>topz
   int_field(k)=int_field(k)+field(k1)*(z_bnd(k1,2)-z_bnd(k1,1));
   z_thick=z_thick+(z_bnd(k1,2)-z_bnd(k1,1));
   if z_bnd(k1,2)==botz
    kix=k1+1;
    break;
   end;
  elseif z_bnd(k1,1)<topz & z_bnd(k1,2)<=botz & z_bnd(k1,2)>topz
   int_field(k)=int_field(k)+field(k1)*(z_bnd(k1,2)-topz); 
   z_thick=z_thick+(z_bnd(k1,2)-topz);
   if  z_bnd(k1,2)==botz
    kix=k1+1;
    break;
   end;
  elseif z_bnd(k1,1)>=topz & z_bnd(k1,2)>botz & z_bnd(k1,1)<=botz
   int_field(k)=int_field(k)+field(k1)*(botz-z_bnd(k1,1)); 
   z_thick=z_thick+(botz-z_bnd(k1,1));
   kix=k1;
   break;
  elseif z_bnd(k1,1)<topz & z_bnd(k1,2)>botz  % if target thickess is inside the layer thickness
   int_field(k)=int_field(k)+field(k1)*(botz-topz);                                
   z_thick=z_thick+(botz-topz);
   kix=k1;
   break;
  elseif z_bnd(k1,1)<topz & z_bnd(k1,2)<topz
   int_field(k)=0.;                                
   kix=k1;
   break;
  end;
  end;
 end;

 if int_field(k)>0.
  int_field(k)=int_field(k)./z_thick;
 end;
end;
% int_field(int_field==0.)=NaN; --> LEO COMMENTED THIS LINE FOR DETOC SINCE SMALL VALUES WERE MADE NaNs in the end of the loop. 
return;


