function [imgr]=dccomp(dctf,q,cf)
[n nn]=size(q);
[M N]=size(dctf);


hb=M/n;
vb=N/n;

dctn=zeros(M,N);
imgr=zeros(M,N);
for i=1:hb
    for j=1:vb
        dctn(n*(i-1)+1:n*(i),n*(j-1)+1:n*j)=dctf(n*(i-1)+1:n*(i),n*(j-1)+1:n*j);
        if(cf==1)
         dctn(n*(i-1)+1:n*(i),n*(j-1)+1:n*j)=dctf(n*(i-1)+1:n*(i),n*(j-1)+1:n*j).*q;
        end
         for x=0:n-1
            for y=0:n-1
                
               
                s=0;
                for u=0:n-1
                    for v=0:n-1
                        au=sqrt(n/2);
                        if(u==0) au=sqrt(n);
                        end
                        av=sqrt(n/2);
                        if(v==0) av=sqrt(n);
                        end
                        au=au*av;
                       t=dctn(n*(i-1)+ u+1,n*(j-1)+ v+1)*cos((pi*(x+0.5)*u)/n)*cos((pi*(y+0.5)*v)/n);
                        t=t/au;
                        s=s+t;
                    end
                end
            
                imgr(n*(i-1)+ x+1,n*(j-1)+ y+1)=s;
            end
        end
       
    end
end
imgr=round(imgr)+128;
imshow(uint8(imgr));
end

