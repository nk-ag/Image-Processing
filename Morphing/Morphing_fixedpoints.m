% running commands-
%prmts
%A2a_Nikita_2018csb1109_2020_CS517('im1.png','im3.png',1,prmts_30,1,'out2wt',1)
function A2a_Nikita_2018csb1109_2020_CS517(fname_inp1,fname_inp2,angle,prmts,afntrms,fname_out1,toshow)
    imgo=imread(fname_inp1);
    imgt=imread(fname_inp2);
    [mo,no]=size(imgo);
    [mt,nt]=size(imgt);
    
    rd=max(ceil((mt-mo)/2),0);
    cd=max(ceil((nt-no)/2),0);
    
    imgo=padarray(imgo,[rd,cd]);
    
    movingpt=[prmts(1,3);prmts(1,4);prmts(2,3);prmts(2,4);prmts(3,3);prmts(3,4);prmts(4,3);prmts(4,4);];
  
    cp=[prmts(:,1)+cd prmts(:,2)+rd];
    mp=[prmts(:,3) prmts(:,4)];
   FileName = sprintf('%s.gif',fname_out1);
    count=1;
   
    for lam=0.0:0.01:1
    
    fname_out=uint8(zeros(mt,nt));
    
   newtie=lam.*cp+(1-lam).*(mp);
     x1=newtie(1,1);y1 =newtie(1,2);
    x2=newtie(2,1);y2 =newtie(2,2);
    x3=newtie(3,1);y3 =newtie(3,2);
    x4=newtie(4,1);y4 =newtie(4,2);
    T=[x1 y1 x1*y1 1 0 0 0 0;
        0 0 0 0 x1 y1 x1*y1 1;
        x2 y2 x2*y2 1 0 0 0 0;
        0 0 0 0 x2 y2 x2*y2 1;
        x3 y3 x3*y3 1 0 0 0 0;
        0 0 0 0 x3 y3 x3*y3 1;
        x4 y4 x4*y4 1 0 0 0 0;
        0 0 0 0 x4 y4 x4*y4 1;];
    c=inv(T)*movingpt;
    for i=1:nt
         for j=1:mt
                 x=round(c(1).*i+c(2).*j+c(3).*(i*j)+c(4));
                 y=round(c(5).*i+c(6).*j+c(7).*(i*j)+c(8));
                
                 if x>0 && y>0 && x<=nt && y<=mt 
                     
                     
                     fname_out(j,i)=imgt(y,x);
                 end
                 
                
         end
          fname_out=fname_out(:,1:nt);

    end
    
   fname_out=trimm(fname_out); %% comment this line to generate image morphing without resize according to transformed image( final image of same size and proportion as original input image)
    if(count==1)
         imwrite(fname_out,FileName,'gif','LoopCount',Inf,'DelayTime',0.05);
        
    else
        imwrite(fname_out,FileName,'gif','WriteMode','append','DelayTime',0.05);
    end
   
    count=count+1;

    end
    
    if(toshow)
        figure;
        subplot(1,3,1);imshow(imgo);title(sprintf('Input image1 %f',lam));subplot(2,3,2);imshow(imgt);title('Input image2');
        subplot(1,3,3);imshow(fname_out);title('Output image');
    end
  

end

function immr=trimm(imrot)
    [d l]=size(imrot);
    f=0;
    for i=1:d
        for j=1:l
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
    for j=1:l
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

    immr=imrot(k:d-k,j:l-j);
    immr=imresize(immr,size(imrot));

end