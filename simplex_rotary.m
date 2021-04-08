%Simplex method code:

%loading variable matrix values
clear
load('TnFile', 'Tn')
load('TdFile', 'Td')
load('ErrorMatrixFile', 'ErrorMatrix')
load('ParametersFile', 'k_d', 'k_u', 'J_2')
%setting parameters for model simulation process
model='two_wheel_rotary';
load_system(model)
set_param(model, 'SimulationMode', 'Normal')

%reflection coefficient
alpha=1;
%compression coefficient
beta=0.5;
%stretching coefficient
gamma=2;

%I. Initial Simplex
Xinit=zeros(1,3);
Yinit=zeros(1,3);
for i=1:3
    Xinit(i)=randi([1100 1300])/100;
    Yinit(i)=randi([1100 1300])/100;
end

%In order to plot the first simplex, we need to know the z-axis values at
%its initial points
errormat=zeros(1,3);
% first initial point error value:
T1=Xinit(1);
T2=Yinit(1);
sim(model)
errormat(1)=error(100001);

% second initial point error value:
T1=Xinit(2);
T2=Yinit(2);
sim(model)
errormat(2)=error(100001);

%third initial point error value:
T1=Xinit(3);
T2=Yinit(3);
sim(model)
errormat(3)=error(100001);

%timer start
tStart=tic;

%surface plot
close all
fig=figure;
fig.Color = [0 0.5 0.5];

surface=surf(Td,Tn,ErrorMatrix, 'FaceAlpha', 0.5);
surface.EdgeColor = 'interp';
colormap(hot)
caxis([0,1.6e8]);

hold on
triangle2=fill3(Yinit,Xinit,errormat, [0 0.8 0.2]);
xlabel('T2')
ylabel('T1')
zlabel('Error Value')

%lowest point is approximately at 820; we know that from scope #2 at the
%end of one code run
plot3(0.02, 2, 820, '.y', 'markersize', 20);

k=0;
LowestError=error(100001);

while LowestError-820>8.2

    %II. Sorting
    
    [errorh, h]=max(errormat);

    X_h=Xinit(h);
    Y_h=Yinit(h);
    Xinit(h)=[];
    Yinit(h)=[];
    errormat(h)=[];
    
    [errorl, l]=min(errormat);
    
    X_l=Xinit(l);
    Y_l=Yinit(l);
    Xinit(l)=[];
    Yinit(l)=[];
    errormat(l)=[];
    
    X_g=Xinit;
    Y_g=Yinit;
    errorg=errormat;
    
    X_center=1/2*sum([X_l,X_g]);
    Y_center=1/2*sum([Y_l, Y_g]);

    X_r=(1+alpha)*X_center-alpha*X_h;
    Y_r=(1+alpha)*Y_center-alpha*Y_h;
    
    T1=X_r;
    T2=Y_r;
    sim(model);
    errorR=error(100001);
    
    %III. Comparison
    
    if errorR<errorl
        X_e=(1-gamma)*X_center+gamma*X_r;
        Y_e=(1-gamma)*Y_center+gamma*Y_r;
        
        T1=X_e;
        T2=Y_e;
        sim(model);
        errore=error(100001);
        
        if errore<errorR
            Xinit=[X_e,X_g,X_l];
            Yinit=[Y_e,Y_g,Y_l];
            errormat(1:3)=[errore, errorg, errorl];
            delete(triangle2)
            triangle2=fill3(Yinit,Xinit,errormat, [0 0.5 0.5]);
            
        elseif errorR<errore 
            Xinit=[X_r,X_g,X_l];
            Yinit=[Y_r,Y_g,Y_l];
            errormat(1:3)=[errorR, errorg, errorl];
            delete(triangle2)
            triangle2=fill3(Yinit,Xinit,errormat, [0.5 0 0.5]);

        end
        
    elseif  errorR<errorg && errorR>errorl
        Xinit=[X_r,X_g,X_l];
        Yinit=[Y_r,Y_g,Y_l];
        errormat(1:3)=[errorR, errorg, errorl];
        delete(triangle2)
        triangle2=fill3(Yinit,Xinit,errormat, [0.5 0.5 0]);

    elseif errorR<errorh && errorR>errorg
        Xinit=[X_r,X_g,X_l];
        Yinit=[Y_r,Y_g,Y_l];
        errormat(1:3)=[errorR, errorg, errorl];
        
        X_s=beta*X_h+(1-beta)*X_center;
        Y_s=beta*Y_h+(1-beta)*Y_center;
        
        T1=X_s;
        T2=Y_s;
        sim(model);
        errorS=error(100001);
        
        if errorS<errorh
            Xinit=[X_s,X_g,X_l];
            Yinit=[Y_s,Y_g,Y_l];
            errormat(1:3)=[errorS, errorg, errorl];
            
            delete(triangle2)
            triangle2=fill3(Yinit,Xinit,errormat, [0.5 0.5 0]);
            
        elseif errorS>=errorh
            
            Xinit=[X_h,X_g,X_l];
            Yinit=[Y_h,Y_g,Y_l];
            errormat(1:3)=[errorh, errorg, errorl];
            
            delete(triangle2)
            triangle2=fill3(Yinit,Xinit,errormat, [0.5 0.5 0]);
            
            for i=1:2
                Xinit(i)=X_l+(Xinit(i)-X_l)/2;
                Yinit(i)=Y_l+(Yinit(i)-Y_l)/2;
                
                T1=Xinit(i);
                T2=Yinit(i);
                sim(model);
                errormat(i)=error(100001);
            end
            delete(triangle2)
            triangle2=fill3(Yinit, Xinit, errormat, [0.25 0.75 0]);
        end
    elseif errorR>errorh
        X_s=beta*X_h+(1-beta)*X_center;
        Y_s=beta*Y_h+(1-beta)*Y_center;
        
        T1=X_s;
        T2=Y_s;
        sim(model);
        errorS=error(100001);
        
        if errorS<errorh
            
            Xinit=[X_s,X_g,X_l];
            Yinit=[Y_s,Y_g,Y_l];
            errormat(1:3)=[errorS, errorg, errorl];
            
            delete(triangle2)
            triangle2=fill3(Yinit, Xinit, errormat, [0.75 0.25 0]);
            
        elseif errorS>errorh
            
            Xinit=[X_h,X_g,X_l];
            Yinit=[Y_h,Y_g,Y_l];
            errormat(1:3)=[errorh, errorg, errorl];
            
            for i=1:2
                Xinit(i)=X_l+(Xinit(i)-X_l)/2;
                Yinit(i)=Y_l+(Yinit(i)-Y_l)/2;
                
                T1=Xinit(i);
                T2=Yinit(i);
                sim(model);
                errormat(i)=error(100001);
               
                
            end
        delete(triangle2)
        triangle2=fill3(Yinit, Xinit, errormat, [0.75 0.25 0]);
           
        end       
    end
    
    k=k+1;
    LowestError=errorl;
end

tSimplex = toc(tStart)
    
