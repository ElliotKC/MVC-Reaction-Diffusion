size = 50;

Ac = zeros(size,size,size); %A concentration
Ac(:) = 1;
Bc = zeros(size,size,size); %B concentration

t = 0; %time
dt = 1; %length of time step

smdmap = linspace(-pi,pi,10);
smd1 = 1;
smd2 = 1;
while smd1 < 11
    while smd2 < 11
    smoothP(smd1,smd2) = (sin(smdmap(smd1)+pi/2)+1)*(sin(smdmap(smd2)+pi/2)+1)/4;
    smd2 = smd2 + 1;
    end
    smd2 = 1;
    smd1 = smd1 + 1;
end

Bc(25:26,25:26,25:26) = 1; %peturbance
%Bc(146:155,146:155) = smoothP;

Da = 1; %A diffusion rate
Db = .5; %b diffusion rate
f = .0545; %feed rate
k = .062; %kill rate

ratioindex = zeros(size,size,floor(t/dt)); %keeps track of data for visualizing

iteration = 0;

while t < 200
    iteration = iteration + 1;
    
    NAc = zeros(size,size,size); %new A concentration
    NAc(:,:) = 1;
    NBc = zeros(size,size,size); %new B concentration

    idx1 = 2;
    idx2 = 2;
    idx3 = 2;

    while idx1 < size
        idx2 = 2;
        while idx2 < size
            idx3 = 2;
            while idx3 < size
                %LaplaceA = Ac(idx1-1,idx2)*.25+Ac(idx1+1,idx2)*.25-Ac(idx1,idx2)+Ac(idx1,idx2-1)*.25+Ac(idx1,idx2+1)*.25;
                LaplaceA = -Ac(idx1,idx2,idx3)+Ac(idx1-1,idx2,idx3)*.2/1.7+Ac(idx1+1,idx2,idx3)*.2/1.7+Ac(idx1,idx2-1,idx3)*.2/1.7+Ac(idx1,idx2+1,idx3)*.2/1.7+Ac(idx1,idx2,idx3-1)*.2/1.7+Ac(idx1,idx2,idx3+1)*.2/1.7+Ac(idx1-1,idx2+1,idx3)*.05/1.7+Ac(idx1-1,idx2-1,idx3)*.05/1.7+Ac(idx1+1,idx2-1,idx3)*.05/1.7+Ac(idx1+1,idx2+1,idx3)*.05/1.7+Ac(idx1-1,idx2,idx3+1)*.05/1.7+Ac(idx1-1,idx2,idx3-1)*.05/1.7+Ac(idx1+1,idx2,idx3-1)*.05/1.7+Ac(idx1+1,idx2,idx3+1)*.05/1.7+Ac(idx1,idx2+1,idx3-1)*.05/1.7+Ac(idx1,idx2-1,idx3-1)*.05/1.7+Ac(idx1,idx2-1,idx3+1)*.05/1.7+Ac(idx1,idx2+1,idx3+1)*.05/1.7+Ac(idx1-1,idx2-1,idx3-1)*0.0125/1.7+Ac(idx1-1,idx2-1,idx3+1)*0.0125/1.7+Ac(idx1-1,idx2+1,idx3-1)*0.0125/1.7+Ac(idx1-1,idx2+1,idx3+1)*0.0125/1.7+Ac(idx1+1,idx2-1,idx3-1)*0.0125/1.7+Ac(idx1+1,idx2-1,idx3+1)*0.0125/1.7+Ac(idx1+1,idx2+1,idx3-1)*0.0125/1.7+Ac(idx1+1,idx2+1,idx3+1)*0.0125/1.7;
                NAc(idx1,idx2,idx3) = Ac(idx1,idx2,idx3)+(Da*LaplaceA-Ac(idx1,idx2,idx3)*Bc(idx1,idx2,idx3)*Bc(idx1,idx2,idx3)+f*(1-Ac(idx1,idx2,idx3)))*dt;

                LaplaceB = -Bc(idx1,idx2,idx3)+Bc(idx1-1,idx2,idx3)*.2/1.7+Bc(idx1+1,idx2,idx3)*.2/1.7+Bc(idx1,idx2-1,idx3)*.2/1.7+Bc(idx1,idx2+1,idx3)*.2/1.7+Bc(idx1,idx2,idx3-1)*.2/1.7+Bc(idx1,idx2,idx3+1)*.2/1.7+Bc(idx1-1,idx2+1,idx3)*.05/1.7+Bc(idx1-1,idx2-1,idx3)*.05/1.7+Bc(idx1+1,idx2-1,idx3)*.05/1.7+Bc(idx1+1,idx2+1,idx3)*.05/1.7+Bc(idx1-1,idx2,idx3+1)*.05/1.7+Bc(idx1-1,idx2,idx3-1)*.05/1.7+Bc(idx1+1,idx2,idx3-1)*.05/1.7+Bc(idx1+1,idx2,idx3+1)*.05/1.7+Bc(idx1,idx2+1,idx3-1)*.05/1.7+Bc(idx1,idx2-1,idx3-1)*.05/1.7+Bc(idx1,idx2-1,idx3+1)*.05/1.7+Bc(idx1,idx2+1,idx3+1)*.05/1.7+Bc(idx1-1,idx2-1,idx3-1)*0.0125/1.7+Bc(idx1-1,idx2-1,idx3+1)*0.0125/1.7+Bc(idx1-1,idx2+1,idx3-1)*0.0125/1.7+Bc(idx1-1,idx2+1,idx3+1)*0.0125/1.7+Bc(idx1+1,idx2-1,idx3-1)*0.0125/1.7+Bc(idx1+1,idx2-1,idx3+1)*0.0125/1.7+Bc(idx1+1,idx2+1,idx3-1)*0.0125/1.7+Bc(idx1+1,idx2+1,idx3+1)*0.0125/1.7;
                NBc(idx1,idx2,idx3) = Bc(idx1,idx2,idx3)+(Da*LaplaceB-Ac(idx1,idx2,idx3)*Bc(idx1,idx2,idx3)*Bc(idx1,idx2,idx3)+f*(1-Bc(idx1,idx2,idx3)))*dt;
                
                if NAc(idx1,idx2,idx3) > 1
                    NAc(idx1,idx2,idx3) = 1;
                end

                if NAc(idx1,idx2,idx3) < 0
                    NAc(idx1,idx2,idx3) = 0;
                end

                if NBc(idx1,idx2,idx3) > 1
                    NBc(idx1,idx2,idx3) = 1;
                end

                if NBc(idx1,idx2,idx3) < 0
                    NBc(idx1,idx2,idx3) = 0;
                end

                idx3 = idx3 + 1;
            end
            idx2 = idx2 + 1;
           
        end
    idx1 = idx1 + 1;
    end
    Ac = NAc;
    Bc = NBc;
    
    figure(1)
    [X,Y,Z] = meshgrid(1:1:50);
    V = Bc./(Ac+Bc);
    xslice = [12.5,25,37.5];   
    yslice = [25];
    zslice = 25;
    slice(X,Y,Z,V,xslice,yslice,zslice)
    view(2)
    shading interp
    t = t+dt;
end

[X,Y,Z] = meshgrid(1:1:50);
V = Bc./(Ac+Bc);
xslice = [12.5,25,37.5];   
yslice = [25];
zslice = 25;
slice(X,Y,Z,V,xslice,yslice,zslice)

