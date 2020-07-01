function [dctf,imgo] =comp(img_p,q)
[n nn]=size(q);
img=imread(img_p);
img=rgb2gray(img);
img=im2double(img);
[M N]=size(img);
rm=mod(M,n);
if(rm~=0) img=[img ;zeros(n-rm,N)];
end
if(rm~=0)
M=M-rm+n;
end
rm=mod(N,n);
if(rm~=0) img=[img zeros(M,n-rm)];
end
if(rm~=0)
N=N-rm+n; end

hb=M/n;
vb=N/n;
imgo=img*255;
img=imgo-128;
r=zeros(n,n);
dctt=zeros(M,N);
for i=1:hb
    for j=1:vb
        for u=0:n-1
            for v=0:n-1
                au=sqrt(n/2);
                if(u==0) au=sqrt(n);
                end
                av=sqrt(n/2);
                if(v==0) av=sqrt(n);
                end
                s=0;
                for x=0:n-1
                    for y=0:n-1
                       s=s+img(n*(i-1)+ x+1,n*(j-1)+ y+1)*cos((pi*(x+0.5)*u)/n)*cos((pi*(y+0.5)*v)/n);
                     
                    end
                end
               au=au*av;
               s=s/au;
                dctt(n*(i-1)+ u+1,n*(j-1)+ v+1)=s;
            end
        end
        dcttf(n*(i-1)+1:n*(i),n*(j-1)+1:n*j)=dctt(n*(i-1)+1:n*(i),n*(j-1)+1:n*j)./q;
    end
end

dctf=round(dcttf);
%imwrite(dctf,"dct90_2.jpg","jpg");
% imwrite(dctt,"dctnr_5.jpg","jpg");
% imwrite(dcr,"dctnl_5.jpg","jpg");

end




