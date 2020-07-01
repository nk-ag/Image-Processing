%Running commands
%prms
%A2b_Nikita_2018csb1109_2020_CS517('manr.jpg','woman.jpeg',tpt,'mg.gif',1)
function A2b_Nikita_2018csb1109_2020_CS517(fname_inp1,fname_inp2,tpts,fname_out,toshow)
img1=imread(fname_inp1);
img2=imread(fname_inp2);
img2=rgb2gray(img2);
FileName=fname_out;

im1=[tpts(:,1),tpts(:,2)];
im2=[tpts(:,3),tpts(:,4)];
%img2=rgb2gray(img2);
 [m1 n1]=size(img1);

%img2=resize(size(img1));
%tri1=delaunayTriangulation(tpts(:,1),tpts(:,2));
tri2=delaunayTriangulation(tpts(:,3),tpts(:,4));

count=1;
for lam=0.0:0.01:1
    dislam=sqrt(1-lam*lam);
    fname_out=uint8(zeros(m1,n1));
    newtie=lam.*im1+(1-lam).*im2;
    tri1=delaunayTriangulation(newtie(:,1),newtie(:,2));
    %tri1.Points()
   % figure;triplot(tri1);
    [ntri, k]=size(tri2.ConnectivityList());
    transf=zeros([6,ntri]);
    
    for i=1:ntri
         pt=tri2.ConnectivityList(i,:);
            mp=tri2.Points(pt,:);
            mp=reshape(mp,[6 1]);
            cp=tri1.Points(pt,:);
            T=zeros([6,6]);
            T(1:3)=cp(1:3);
            T(22:24)=cp(1:3);
            T(7:9)=cp(4:6);
            T(28:30)=cp(4:6);
            T(13:15)=uint8(ones([3 1]));
            T(34:36)=uint8(ones([3 1]));
            c=inv(T)*mp;
            transf(:,i)=c;
            
    end
    for i=1:m1
        for j=1:n1
            id=pointLocation(tri2,[i,j]);
            c=transf(:,id);
            x=round(c(1).*i+c(2).*j+c(3));
             y=round(c(4).*i+c(5).*j+c(6));

             if x>0 && y>0 && x<=m1 && y<=n1 


                 fname_out(j,i)=(dislam)*img1(y,x)+(1-dislam)*(img2(y,x));
             end
        end
        %fname_out=fname_out(:,1:n1);
    end
    
    if(count==1)
         imwrite(fname_out,FileName,'gif','LoopCount',Inf,'DelayTime',0.05);
        
    else
        imwrite(fname_out,FileName,'gif','WriteMode','append','DelayTime',0.05);
    end
   
    count=count+1;

end
 imwrite(fname_out,FileName,'gif','WriteMode','append','DelayTime',0.5);
    
     if(toshow)
        figure;
        subplot(1,4,1);imshow(img1);title(sprintf('Input image1 %f',lam));subplot(2,4,2);imshow(img2);title('Input image2');
        subplot(1,4,3);imshow(fname_out);title('Output image');subplot(1,4,4);imshow(fname_out-img2);title('Error image');
    end


