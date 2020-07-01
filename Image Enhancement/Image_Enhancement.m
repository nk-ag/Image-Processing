function outx=A1_Nikita_2018csb1109_2019_CS517(qID,fname_inp1,fname_inp2,fname_out,prmts,toshow)
    outx=[];
    if(qID==1), outx=interp_NN(fname_inp1,fname_out,prmts,toshow);end
    if(qID==2), outx=interp_bilinear(fname_inp1,fname_out,prmts,toshow);end
    if(qID==3), outx=rotateim(fname_inp1,fname_out,prmts,toshow);end
    if(qID==4), outx=bitplane(fname_inp1,fname_out,prmts,toshow);end
    if(qID==5), outx=tie_points(fname_inp1,fname_inp2,fname_out,prmts,toshow);end
    if(qID==6), outx=histo_eq(fname_inp1,fname_out,toshow);end
    if(qID==7), outx=histo_match(fname_inp1,fname_inp2,fname_out,toshow);end
    if(qID==8), outx=histo_adaptive(fname_inp1,fname_out,toshow);end
    
end

function rmse=interp_bilinear(fname_inp1,fname_out,prmts,toshow)
img=imread(fname_inp1);
m1=prmts(1);
n1=prmts(2);
[m,n]=size(img);
rsc=m/m1;
csc=n/n1;
fname_out=uint8(zeros(m1,n1));
SImg=imresize(img,[m1,n1],'bilinear');

for i=1:m1
   for j=1:n1
        i1=(i).*rsc;
        j1=(j).*csc;
        i1f=max(floor(i1),1);
        j1f=max(floor(j1),1);
        i1c=min(i1f+1,m);
        j1c=min(j1f+1,n);
        p1=img(i1f,j1f).*(i1-i1f)+(img(i1c,j1f).*(i1c-i1));
        p2=img(i1f,j1c).*(i1-i1f)+(img(i1c,j1c).*(i1c-i1));
        p3=p1.*max(j1c-j1)+p2.*(j1-j1f);
        fname_out(i,j)=p3;
   end 
end

rmse=sqrt(immse(SImg,fname_out));
if(toshow)
    figure;
    subplot(2,2,1);imshow(img);title('Input Image');subplot(2,2,2);imshow(fname_out);title('Output Image');subplot(2,2,3);imshow(SImg);title('SystemImage');subplot(2,2,4);imshow(imresize(img,[m1,n1],'bicubic'));title('Bicubic');
end

end

function rmse=interp_NN(fname_inp1,fname_out,prmts,toshow)
    img=imread(fname_inp1);
    m1=prmts(1);
    n1=prmts(2);
    [m,n]=size(img);
    rsc=m/m1;
    csc=n/n1;
    fname_out=uint8(zeros(m1,n1));
    SImg=imresize(img,[m1,n1],'nearest');
    for i=1:m1
       for j=1:n1
            i1=round(0.5+((i-0.5)*rsc));
            j1=round(0.5+((j-0.5)*csc));
            fname_out(i,j)=img(i1,j1);
       end
    end
    rmse=sqrt(immse(SImg,fname_out));
    if(toshow)
        figure;
        subplot(2,2,1);imshow(img);title('Input Image');subplot(2,2,2);imshow(fname_out);title('Output Image');subplot(2,2,3);imshow(SImg);title('SystemImage');subplot(2,2,4);imshow(imresize(img,[m1,n1],'bicubic'));title('Bicubic');
    end

end

function o_dim=rotateim(fname_inp1,fname_out,prmts,toshow)
    img=imread(fname_inp1);
    th=prmts;
    T=[cosd(th),sind(th),0; -sind(th),cosd(th),0;0,0,1];
    [m,n]=size(img);
    d=ceil(sqrt(m.^2+n.^2));
    imrot=uint8(zeros(d,d));
    imgp=uint8(zeros(d,d));
    i=ceil((d+m)/2)-floor((d-m)/2)-m+1;
    j=ceil((d+n)/2)-floor((d-n)/2)-n+1;
    imgp(floor((d-m)/2):ceil((d+m)/2)-i,floor((d-n)/2):ceil((d+n)/2)-j)=img;
    midx=ceil((size(imgp,1)+1)/2);
    midy=ceil((size(imgp,2)+1)/2);
    for i=1:d
        for j=1:d
            k=[i-midx j-midy 1]/T;
            x=round(k(1))+midx;
            y=round(k(2))+midy;
            if x>0 && y>0 && x<=d &&y<=d
                imrot(i,j)=imgp(x,y);
            end
        end
    end
    f=0;
    for i=1:d
        for j=1:d
            if(imrot(i,j)~=0)
                f=1;
                break;
            end
        end
       if(f==1)
           break;
       end
    end

    k=i;
    f=0;
    for j=1:d
        for i=1:d
            if(imrot(i,j)~=0)
                f=1;
                break;
            end
        end
       if(f==1)
           break;
       end
    end

    fname_out=imrot(k:d-k,j:d-j);
    SImg=imrotate(img,th);
    t=imresize(SImg,size(fname_out));
    rmse=sqrt(immse(fname_out,t));
    o_dim=[size(fname_out),rmse];
    if(toshow)
        figure;
        subplot(1,3,1);imshow(img);title('Input Image');subplot(1,3,2);imshow(fname_out);title('Output Image');subplot(1,3,3);imshow(SImg);title('System Image');
    end
    


end

function rmse=bitplane(fname_inp1,fname_out,prmts,toshow)
img=imread(fname_inp1);
[m, n]=size(img);
bpx=prmts;
[~,l]=size(prmts);
figure;
rmse=[];
for i=1:l
    b=bpx(i);
    bin=tobinary(b);
    fname_out=uint8(zeros(m,n));
    for p=1:m
        for j=1:n
            t=tobinary(img(p,j));
            k=255-todecimal(t&bin);
            fname_out(p,j)=k;
        end
    end
    rmse(i)=sqrt(immse(img,fname_out));
    
    
    if(toshow)
        
        subplot(1,l+1,i+1);imshow(fname_out);title('Output image');
    end
end
subplot(1,l,1);imshow(img);


end

function bin= tobinary(num)
    bin=[0 0 0 0 0 0 0 0 ];
    for i=1:8
        
        bin(i)=mod(num,2);
        num=floor(num/2);
        
    end
    
end

function dec=todecimal(bin)
dec=0;

    for i=0:7
        dec=2.^i*bin(i+1)+dec;
        
    end
    

end


function o_dim=tie_points(fname_inp1,fname_inp2,fname_out,prmts,toshow)
    imgo=imread(fname_inp1);
    imgt=imread(fname_inp2);
    [mo,no]=size(imgo);
    [mt,nt]=size(imgt);
    movingpt=[prmts(1,3);prmts(1,4);prmts(2,3);prmts(2,4);prmts(3,3);prmts(3,4);prmts(4,3);prmts(4,4);];
    x1=prmts(1,1);y1 =prmts(1,2);
    x2=prmts(2,1);y2 =prmts(2,2);
    x3=prmts(3,1);y3 =prmts(3,2);
    x4=prmts(4,1);y4 =prmts(4,2);
    T=[x1 y1 x1*y1 1 0 0 0 0;
        0 0 0 0 x1 y1 x1*y1 1;
        x2 y2 x2*y2 1 0 0 0 0;
        0 0 0 0 x2 y2 x2*y2 1;
        x3 y3 x3*y3 1 0 0 0 0;
        0 0 0 0 x3 y3 x3*y3 1;
        x4 y4 x4*y4 1 0 0 0 0;
        0 0 0 0 x4 y4 x4*y4 1;];
    
    c=inv(T)*movingpt;
    fname_out=uint8(zeros(mo,no));
    for i=1:mo
         for j=1:no
                 x=round(c(1).*i+c(2).*j+c(3).*(i*j)+c(4));
                 y=round(c(5).*i+c(6).*j+c(7).*(i*j)+c(8));
                 if x>0 && y>0 && x<=nt && y<=mt
                     %print(i,j);
                     fname_out(j,i)=imgt(y,x);
                 end
         end
    end
    rmse=sqrt(immse(imgo,fname_out));
    o_dim=[size(fname_out),rmse];
    %figure;imshow(fname_out);
    if(toshow)
        figure;
        subplot(2,2,1);imshow(imgo);title('Input image1');subplot(2,2,2);imshow(imgt);title('Input image2');
        subplot(2,2,3);imshow(fname_out);title('Output image');subplot(2,2,4);imshow(imgo-fname_out);title('Error Image');
    end

end

function rmse= histo_eq(fname_inp1,fname_out,toshow)

    I=imread(fname_inp1);
    [m,n]=size(I);
    freq=zeros(256,1);
    cfreq=zeros(256,1);
    out=zeros(256,1);
    fname_out=uint8(zeros(m,n));
    for i=1:m
        for j=1:n
            freq(I(i,j)+1)=freq(I(i,j)+1)+1;

        end
    end
    cfreq(1)=freq(1);
    for i=1:255
        cfreq(i+1)=cfreq(i)+freq(i+1);
        out(i+1)=ceil(255.*(cfreq(i+1)));
    end
    out=out./cfreq(end);
    for i=1:m
        for j=1:n
            fname_out(i,j)=out(I(i,j)+1);

        end

    end
    SImg=histeq(I);
    rmse=sqrt(immse(SImg,fname_out));
    if(toshow)
        figure;
        subplot(2,3,1);imshow(I);title('Input Image');subplot(2,3,2);imshow(fname_out);title('Output Image');subplot(2,3,3);imshow(SImg);title('System Image');
        subplot(2,3,4);imhist(I);title('Input Image Histogram');subplot(2,3,5);imhist(fname_out);title('Output Image histogram');subplot(2,3,6);imhist(SImg);title('System Image Histogram');

    end
end

function rmse=histo_match(fname_inp1,fname_inp2,fname_out,toshow)
imr=imread(fname_inp1);
imz=imread(fname_inp2);
freqr=zeros(256,1);
cfreqr=zeros(256,1);
freqz=zeros(256,1);
cfreqz=zeros(256,1);
[mr ,nr]=size(imr);
[mz ,nz]=size(imz);
for i=1:mr
    for j=1:nr
        freqr(imr(i,j)+1)=freqr(imr(i,j)+1)+1;
    end
end
for i=1:mz
    for j=1:nz
        freqz(imz(i,j)+1)=freqz(imz(i,j)+1)+1;
    end
end
cfreqr(1)=freqr(1);
cfreqz(1)=freqz(1);
for i=2:256
    cfreqr(i)=cfreqr(i-1)+freqr(i);
    cfreqz(i)=cfreqz(i-1)+freqz(i);
    
end
cfreqr=cfreqr./cfreqr(end);
cfreqz=cfreqz./cfreqz(end);
cfreqr=255.*cfreqr;
cfreqz=255.*cfreqz;
inv=zeros(256,1);

for i=1:256
    s=cfreqr(i);
    f=0;
   % inv(i)=cfreqz(1);
    for j=1:255
        if cfreqz(j+1)>s && cfreqz(j)<s 
            f=1;
            if(cfreqz(j+1)==s)
                inv(i)=j;
            elseif(cfreqz(j)==s)
                    inv(i)=j-1;
            else
                if (cfreqr(j+1)-s)<(s-cfreqr(j))
                     inv(i)=j;
                
                else
                    inv(i)=j-1;
                end
            end
            break;
        end
        
    end
        if(f==0)
            inv(i)=255;
        end
end
fname_out=uint8(zeros(mr, nr));
for i=1:mr
    for j=1:nr
        fname_out(i, j)=inv(imr(i,j)+1);
    
    end

end

if(toshow)
    figure;
    subplot(2,3,1);imshow(imr);title('Input Image 1');subplot(2,3,2);imshow(imz);title('Input Image 2');subplot(2,3,3);imshow(fname_out);title('Output Image');
    subplot(2,3,4);imhist(imr);title('Input Image 1 Histogram');subplot(2,3,5);imhist(imz);title('Input Image 2 Histogram');subplot(2,3,6);imhist(fname_out);title('Output Image Histogram');

end
end

function rmse=histo_adaptive(fname_inp1,fname_out,toshow)
img=imread(fname_inp1);
[m, n]=size(img);
fname_out=uint8(zeros(m,n));
for i=2:m-2
    
    

    for j=2:n-2
        freq=zeros(256,1);
        cfreq=zeros(256,1);
        for k=-1:1
            for t=-1:1
            freq(img(i+k,j+t)+1)=freq(img(i+k,j+t)+1)+1;
            end
        end
        cfreq(1)=freq(1);
        for k=2:256
            cfreq(k)=cfreq(k-1)+freq(k);
            
        end
        cfreq=cfreq./cfreq(end);
        cfreq=ceil(255.*cfreq);
        fname_out(i,j)=cfreq(img(i,j)+1);

    end

end
sysHeq=adapthisteq(img);
rmse=sqrt(immse(sysHeq,fname_out));
if(toshow)
    figure;
    subplot(2,3,1);imshow(img);title('Input Image');subplot(2,3,2);imshow(fname_out);title('Output Image');subplot(2,3,3);imshow(sysHeq);title('System Image');
    subplot(2,3,4);imhist(img);title('Input Image Histogram');subplot(2,3,5);imhist(fname_out);title('Output Image Histogram');subplot(2,3,6);imhist(sysHeq);title('Systems Image Histogram');
end


end