clear 
clc 

sigma = [0.0018 0.0027 0.0008;
         0.0027 0.0102 0.0048;
         0.0008 0.0048 0.0433]; 
     
R = [0.006 0.021 0.09]'; 

V = @(w) w'*sigma*w; 

%remember that the default for most of matlab functions use the convention
%That vectors are column vectors 

value = fmincon(V,[1 1 1]',[],[],[1 1 1;R'],[1;0.07],[0 0 0]',[0.1 0.2 1]);
